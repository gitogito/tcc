#!/usr/bin/env ruby

require 'eggx'

class PPM
  attr_reader :width, :height
  attr_reader :rgb_ary

  def initialize(io)
    @io = io
    @tokens = []

    unless next_token == 'P3'
      raise 'invalid ppm file'
    end
    @width = Integer(next_token)
    @height = Integer(next_token)
    depth = Integer(next_token)
    raise 'depth != 255' unless depth == 255
    @rgb_ary = Array.new
    loop do
      token1 = next_token
      break if token1.nil?
      token2 = next_token
      token3 = next_token
      @rgb_ary <<
        (Integer(token1) << 16) + (Integer(token2) << 8) + Integer(token3)
    end
  end

  private

  def next_token
    while @tokens.empty?
      line = @io.gets
      next if /^#/ =~ line
      if line.nil?
        return nil
      end
      @tokens = line.split
    end
    @tokens.shift
  end
end

def redraw(win, ppm, xy_ary)
  win.putimg24(0.0, 0.0, ppm.width, ppm.height, ppm.rgb_ary)
  unless xy_ary.empty?
    win.line(xy_ary[0][0], xy_ary[0][1], EGGX::PENUP)
    xy_ary.each do |x, y|
      win.line(x, y, EGGX::PENDOWN)
    end
  end
end

ppm = PPM.new(ARGF)
win = EGGX::Window.new(ppm.width, ppm.height)
win.putimg24(0.0, 0.0, ppm.width, ppm.height, ppm.rgb_ary)

xy_ary = []
x, y = 0, 0
manhattan = false
colors = [[255, 255, 0], [0, 255, 255], [255, 0, 255], [255, 255, 255], [0, 0, 0]]
win.newrgbcolor(*colors[0])
loop do
  inp = win.getxpress
  case inp[0]
  when EGGX::ButtonPress
    x, y = inp[2, 2]
    if manhattan and not xy_ary.empty?
      if (x - xy_ary[-1][0]).abs > (y - xy_ary[-1][1]).abs
        y = xy_ary[-1][1]
      else
        x = xy_ary[-1][0]
      end
    end
    if xy_ary.empty?
      win.line(x, y, EGGX::PENUP)
    else
      win.line(x, y, EGGX::PENDOWN)
    end
    xy_ary << [x, y]
    STDERR.print '(', x, ', ', y, ')', "\n"
  when EGGX::KeyPress
    case inp[1]
    when ?c
      colors = colors[1 .. -1] + [ colors[0] ]
      win.newrgbcolor(*colors[0])
      redraw(win, ppm, xy_ary)
    when ?d
      xy_ary.pop
      redraw(win, ppm, xy_ary)
    when ?m
      manhattan = ! manhattan
      STDERR.puts "manhattan is #{manhattan}"
    when ?q
      break
    end
  end
end

xy_ary.each do |x, y|
  print x, ', ', y, "\n"
end
