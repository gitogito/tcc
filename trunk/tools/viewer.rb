#!/usr/bin/env ruby

require 'narray'
require 'eggx'
require 'optparse'

WINNAME = 'viewer'
BLOCK = 10
#COLOR = EGGX::IDL2_BLUE_RED
#COLOR = EGGX::DS9_GRAY
COLOR = EGGX::DS9_B
#COLOR = EGGX::IDL1_RED_TEMPERATURE
FONTMARGIN = 2
FONTSIZE = 16
BGCOLOR = 'white'
TEXTCOLOR = 'black'
BORDERCOLOR = 'blue'

opts = Hash.new
OptionParser.new do |parser|
  parser.on('-o outfile') do |arg|
    opts['o'] = String(arg)
  end

  begin
    parser.parse!
  rescue => e
    STDERR.puts e.message
    STDERR.puts parser.help
    exit 1
  end
end

outfile = (opts['o'] or nil)

if gets =~ /\A#\s*(\S+)\s+(\S+)\s+(\S+)\s*\Z/
  nx, xmin, xmax = [Integer($1), Float($2), Float($3)]
else
  raise "can't read"
end
if gets =~ /\A#\s*(\S+)\s+(\S+)\s+(\S+)\s*\Z/
  ny, ymin, ymax = [Integer($1), Float($2), Float($3)]
else
  raise "can't read"
end
if gets =~ /\A#\s*(\S+)\s+(\S+)\s+(\S+)\s*\Z/
  nz, zmin, zmax = [Integer($1), Float($2), Float($3)]
else
  raise "can't read"
end

x3d = NArray.float(nx, ny, nz)
y3d = NArray.float(nx, ny, nz)
z3d = NArray.float(nx, ny, nz)
u3d = NArray.float(nx, ny, nz)
umin = umax = nil
nz.times do |k|
  ny.times do |j|
    nx.times do |i|
      line = gets
      redo if line =~ /\A\Z/ or line =~ /\A#/
      x, y, z, u = line.split.map{|v| Float(v)}
      x3d[i, j, k] = x
      y3d[i, j, k] = y
      z3d[i, j, k] = z
      u3d[i, j, k] = u

      umin = u if umin.nil? or umin > u
      umax = u if umax.nil? or umax < u
    end
  end
end

step = 1
block = BLOCK
draw_p = true
win = nil
xwin = ywin = 0
surface = :ZX
i_a3 = 0
draw_width = draw_height = nil
border = false
loop do
  case surface
  when :YZ
    a1_2d = y3d[i_a3, nil, nil]
    a2_2d = z3d[i_a3, nil, nil]
    u2d   = u3d[i_a3, nil, nil]
    a3    = x3d[i_a3, 0  , 0  ]
    a1min = ymin
    a1max = ymax
    a2min = zmin
    a2max = zmax
    n_a3 = nx
  when :ZX
    a1_2d = x3d[nil, i_a3, nil]
    a2_2d = z3d[nil, i_a3, nil]
    u2d   = u3d[nil, i_a3, nil]
    a3    = y3d[0  , i_a3, 0  ]
    a1min = xmin
    a1max = xmax
    a2min = zmin
    a2max = zmax
    n_a3 = ny
  when :XY
    a1_2d = x3d[nil, nil, i_a3]
    a2_2d = y3d[nil, nil, i_a3]
    u2d   = u3d[nil, nil, i_a3]
    a3    = z3d[0  , 0  , i_a3]
    a1min = xmin
    a1max = xmax
    a2min = ymin
    a2max = ymax
    n_a3 = nz
  else
    raise 'bug'
  end

  n1, n2 = u2d.shape
  if draw_p
    draw_width = block * n1
    draw_height = block * n2
    EGGX.setinitialwinname(WINNAME, WINNAME, WINNAME, WINNAME)
    EGGX.setinitialbgcolor(BGCOLOR)
    win.close unless win.nil?
    win = EGGX::Window.new(draw_width, draw_height + 2*(FONTMARGIN + FONTSIZE))
    draw_p = false
  end
  win.layer(0, 1)
  if border
    win.newcolor(BORDERCOLOR)
    win.fillrect(0, 0, draw_width, draw_height)
  end
  n2.times do |j|
    n1.times do |i|
      r, g, b = EGGX.makecolor(COLOR, umin, umax, u2d[i, j])
      win.newrgbcolor(r, g, b) ;
      if border
	win.fillrect(block * i, block * j, block-1, block-1)
      else
	win.fillrect(block * i, block * j, block, block)
      end
    end
  end
  x = a1min + xwin.to_f / draw_width * (a1max - a1min)
  y = a2min + ywin.to_f / draw_height * (a2max - a2min)
  i = ((x - a1min) / (a1max - a1min) * n1).to_i
  j = ((y - a2min) / (a2max - a2min) * n2).to_i
  win.newcolor(BGCOLOR)
  win.fillrect(0, draw_height, block * n1, 2*(FONTMARGIN+FONTSIZE))
  win.newcolor(TEXTCOLOR)
  win.drawstr(0, block * n2 + FONTMARGIN, FONTSIZE, 0,
  "%.3e: %.1e, %.1e" % [u2d[i, j], a1_2d[i, j], a2_2d[i, j]])
  win.drawstr(0, block * n2 + FONTMARGIN+FONTSIZE+FONTMARGIN, FONTSIZE, 0,
  "%s: %.2e: %.1f%" % [surface.to_s, a3, i_a3.to_f / (n_a3 - 1) * 100])
  win.copylayer(1, 0)

  type, code, xwin_tmp, ywin_tmp = win.getxpress
  case type
  when EGGX::ButtonPress
    if xwin_tmp >= 0 and xwin_tmp < draw_width and
      ywin_tmp >= 0 and ywin_tmp < draw_height
      xwin = xwin_tmp
      ywin = ywin_tmp
    end
  when EGGX::KeyPress
    case code
    when ?q.ord
      break
    when ?s.ord
      if outfile.kind_of?(String)
	f = open(outfile, 'w')
      else
	f = STDOUT
      end
      f.puts "# %s: %.2e: %.1f%" % [surface.to_s, a3, i_a3.to_f / (n_a3 - 1) * 100]
      n2.times do |j|
	n1.times do |i|
	  x = a1min + (a1max - a1min) / (n1-1) * i
	  y = a2min + (a2max - a2min) / (n2-1) * j
	  f.puts "%g\t%g\t%g" % [x, y, u2d[i, j]]
	end
	f.print "\n"
      end
      f.close unless f == STDOUT
      break
    when ?x.ord
      surface = :YZ
      i_a3 = 0
      xwin = 0
      ywin = 0
      draw_p = true
    when ?y.ord
      surface = :ZX
      i_a3 = 0
      xwin = 0
      ywin = 0
      draw_p = true
    when ?z.ord
      surface = :XY
      i_a3 = 0
      xwin = 0
      ywin = 0
      draw_p = true
    when ?j.ord
      if i_a3 >= step
	i_a3 -= step
      else
	i_a3 = 0
      end
    when ?k.ord
      if i_a3 < n_a3 - step
	i_a3 += step
      else
	i_a3 = n_a3 - 1
      end
    when ?1.ord; step = 1
    when ?2.ord; step = 2
    when ?3.ord; step = 3
    when ?4.ord; step = 4
    when ?5.ord; step = 5
    when ?>.ord
      block += step
      draw_p = true
      xwin = 0
      ywin = 0
    when ?<.ord
      if block - step >= 1
	block -= step
	draw_p = true
      else
        block = 1
	draw_p = true
      end
      xwin = 0
      ywin = 0
    when ?b.ord
      if border
	border = false
      else
	border = true
      end
    end
  end
end
