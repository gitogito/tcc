module TCC
  class Model
    def initialize(pnt, vec, d)
      puts('--- world')
      puts('(%g, %g, %g), <%g, %g, %g>, %g, %g, %g' % [
           pnt[0], pnt[1], pnt[2],
           vec[0], vec[1], vec[2],
           d[0], d[1], d[2]])
      puts('')
    end

    def active(obj)
      puts('--- active')
      puts obj.to_s
      puts('')
    end

    def inactive(obj)
      puts('--- inactive')
      puts obj.to_s
      puts('')
    end

    def fix(temp, obj)
      puts('--- fix %g' % [temp])
      puts obj.to_s
      puts('')
    end

    def fixheat(power, obj)
      puts('--- fix heat %g' % [power])
      puts obj.to_s
      puts('')
    end

    def heat(power, obj)
      puts('--- heat %g' % [power])
      puts obj.to_s
      puts('')
    end

    def lambda(k, obj)
      puts('--- lambda %g' % [k])
      puts obj.to_s
      puts('')
    end
  end

  class Obj
    IND = '  '
    @@nind = 0

    private

    def indent
      IND * @@nind
    end
  end

  class Line < Obj
    def initialize(axis, pnt, vec2d_ary)
      @axis = axis
      @pnt = pnt
      @vec2d_ary = vec2d_ary
    end

    def to_s
      s = indent
      s += 'line :%s, (%g, %g, %g)' % [@axis, @pnt[0], @pnt[1], @pnt[2]]
      @vec2d_ary.each do |vec2d|
        s += ', <%g, %g>' % [vec2d[0], vec2d[1]]
      end
      s
    end
  end

  class Polygon < Obj
    def initialize(axis, pnt, vec2d_ary)
      @axis = axis
      @pnt = pnt
      @vec2d_ary = vec2d_ary
    end

    def to_s
      s = indent
      s += 'polygon :%s, (%g, %g, %g)' % [@axis, @pnt[0], @pnt[1], @pnt[2]]
      @vec2d_ary.each do |vec2d|
        s += ', <%g, %g>' % [vec2d[0], vec2d[1]]
      end
      s
    end
  end

  class Triangle < Obj
    def initialize(axis, pnt, vec2d_a, vec2d_b)
      @axis = axis
      @pnt = pnt
      @vec2d_a = vec2d_a
      @vec2d_b = vec2d_b
    end

    def to_s
      s = indent
      s += sprintf('triangle :%s, (%g, %g, %g), <%g, %g>, <%g, %g>',
                   @axis,
                   @pnt[0], @pnt[1], @pnt[2],
                   @vec2d_a[0], @vec2d_a[1],
                   @vec2d_b[0], @vec2d_b[1])
      s
    end
  end

  class Rect < Obj
    def initialize(axis, pnt, vec2d)
      @axis = axis
      @pnt = pnt
      @vec2d = vec2d
    end

    def to_s
      s = indent
      s += sprintf('rect :%s, (%g, %g, %g), <%g, %g>',
                   @axis,
                   @pnt[0], @pnt[1], @pnt[2],
                   @vec2d[0], @vec2d[1])
      s
    end
  end

  class EllipsePeri < Obj
    NDIV = 50

    def initialize(axis, pnt, ru, rv, st_angle, en_angle, polygon_p = false)
      st_theta = 2 * Math::PI / 360.0 * st_angle
      en_theta = 2 * Math::PI / 360.0 * en_angle
      u0 = ru * Math.cos(st_theta)
      v0 = rv * Math.sin(st_theta)
      case axis
      when :X
        pnt0 = [pnt[0]     , pnt[1] + u0, pnt[2] + v0]
      when :Y
        pnt0 = [pnt[0] + u0, pnt[1]     , pnt[2] + v0]
      when :Z
        pnt0 = [pnt[0] + u0, pnt[1] + v0, pnt[2]     ]
      else
        raise "unknown axis #{axis}"
      end
      vec2d_ary = []
      (1 ... NDIV).each do |i|
        theta = st_theta + Float(en_theta - st_theta) / NDIV * i
        u = ru * Math.cos(theta)
        v = rv * Math.sin(theta)
        vec2d_ary << [u - u0, v - v0]
      end
      if polygon_p
        @obj = Polygon.new(axis, pnt0, vec2d_ary)
      else
        @obj = Line.new(axis, pnt0, vec2d_ary)
      end
    end

    def to_s
      @obj.to_s
    end
  end

  class CirclePeri < Obj
    def initialize(axis, pnt, radius, st_angle, en_angle)
      @ellipse_peri = EllipsePeri.new(axis, pnt, radius, radius,
                                      st_angle, en_angle)
    end

    def to_s
      @ellipse_peri.to_s
    end
  end

  class Ellipse < Obj
    def initialize(axis, pnt, ru, rv)
      @ellipse_peri = EllipsePeri.new(axis, pnt, ru, rv, 0.0, 360.0, true)
    end

    def to_s
      @ellipse_peri.to_s
    end
  end

  class Circle < Obj
    def initialize(axis, pnt, radius)
      @ellipse = Ellipse.new(axis, pnt, radius, radius)
    end

    def to_s
      @ellipse.to_s
    end
  end

  class Box < Obj
    def initialize(pnt, vec)
      @pnt = pnt
      @vec = vec
    end

    def to_s
      s = indent
      s += sprintf('box (%g, %g, %g), <%g, %g, %g>',
                   @pnt[0], @pnt[1], @pnt[2],
                   @vec[0], @vec[1], @vec[2])
      s
    end
  end

  class Objs < Obj
    def initialize(objs)
      @objs = objs
    end

    def self.[](*objs)
      self.new(objs)
    end

    def <<(obj)
      @objs << obj
    end

    def to_s
      s = indent
      s += "[\n"
      @@nind += 1
      (@objs.size - 1).times do |i|
        s += @objs[i].to_s
        s += "\n"
      end
      s += @objs[-1].to_s
      s += "\n"
      @@nind -= 1
      s += indent
      s += "]"
      s
    end
  end

  class Sweep < Obj
    def initialize(axis, length, obj)
      @axis = axis
      @length = length
      @obj = obj
    end

    def to_s
      s = indent
      s += sprintf("sweep :%s, %g,\n", @axis, @length)
      @@nind += 1
      s += @obj.to_s
      @@nind -= 1
      s
    end
  end
end

if __FILE__ == $0
  include TCC

  tcc = Model.new([0, 0, 0], [1, 1, 1], [1.0/30, 1.0/30, 1.0/30])
  tcc.active(
    Objs[
      Sweep.new(:Z, 0.5,
                Rect.new(:Z, [0, 0, 0], [1, 1])),
      Box.new([0, 0, 0.5], [0.5, 0.5, 0.5])
  ])
  tcc.fix(0, Rect.new(:Z, [0, 0, 0], [1, 1]))
  tcc.fix(1, Rect.new(:Z, [0, 0, 1], [0.5, 0.5]))
  tcc.heat(1, Circle.new(:Z, [0.5, 0.5, 0.2], 0.2))
  tcc.fix(0.5, Ellipse.new(:Z, [0.5, 0.5, 0.1], 0.2, 0.3))
  tcc.lambda(0.1, Box.new([0, 0, 0.5], [0.5, 0.5, 0.5]))
  tcc.fix(0.9, Line.new(:Z, [0, 0, 0.3], [[0.5, 0.5], [0.1, 0.8]]))
  tcc.inactive(Sweep.new(:Z, 0.5,
                         Circle.new(:Z, [0.5, 0.5, 0.5], 0.2)))
  tcc.active(Sweep.new(:Z, 0.2,
                       Polygon.new(:Z, [0, 0, 0.8], [[1, 0], [1, 1]])))
  tcc.active(Sweep.new(:Z, 0.2,
                       Triangle.new(:Z, [0.5, 0.5, 0.8], [0.1, 0.5], [-0.2, 0.5])))
end
