class TcParser

options no_result_var

prechigh
  right TK_POW
  right NEG
  left '*' '/'
  left '+' '-'
preclow

rule

input:
    world commands

  | var_assigns world commands

var_assigns:
    var_assign

  | var_assigns var_assign

var_assign:
    TK_WORD '=' expr
    { @vars[val[0]] = val[2] }

world:
    TK_DASH TK_WORLD point ',' vector ',' expr ',' expr ',' expr
    {
      obj = Box.new(*val.values_at(2, 4))
      obj.type = :WORLD
      @obj_ary << obj
    }

  | TK_DASH TK_WORLD point ',' point ',' expr ',' expr ',' expr
    {
      val[4][0] -= val[2][0]
      val[4][1] -= val[2][1]
      val[4][2] -= val[2][2]
      obj = Box.new(*val.values_at(2, 4))
      obj.type = :WORLD
      @obj_ary << obj
    }

commands:
    command

  | commands command

command:
    TK_DASH TK_ACTIVE obj
    {
      val[2].type = :ACTIVE
      @obj_ary << val[2]
    }

  | TK_DASH TK_INACTIVE obj
    {
      val[2].type = :INACTIVE
      @obj_ary << val[2]
    }

  | TK_DASH TK_FIX expr obj
    {
      val[3].type = :FIX
      @obj_ary << val[3]
    }

  | TK_DASH TK_FIX TK_HEAT expr obj
    {
      val[4].type = :FIXHEAT
      @obj_ary << val[4]
    }

  | TK_DASH TK_HEAT TK_FIX expr obj
    {
      val[4].type = :FIXHEAT
      @obj_ary << val[4]
    }

  | TK_DASH TK_HEAT expr obj
    {
      val[3].type = :HEAT
      @obj_ary << val[3]
    }

  | TK_DASH TK_LAMBDA expr obj
    {
      val[3].type = :LAMBDA
      @obj_ary << val[3]
    }

point:
    '(' expr ',' expr ',' expr ')'
    {
      point = val.values_at(1, 3, 5)
      set_min_max(point)
      point
    }

point2d:
    '(' expr ',' expr ')'
    {
      point2d = val.values_at(1, 3)
      set_min_max(point2d)
      point2d
    }

point2d_ary:
    point2d
    { [ val[0] ] }

  | point2d_ary ',' point2d
    { val[0] + [ val[2] ] }

vector:
    '<' expr ',' expr ',' expr '>'
    {
      vector = val.values_at(1, 3, 5)
      set_min_max(vector)
      vector
    }

vector2d:
    '<' expr ',' expr '>'
    {
      vector2d = val.values_at(1, 3)
      set_min_max(vector2d)
      vector2d
    }

vector2d_ary:
    vector2d
    { [ val[0] ] }

  | vector2d_ary ',' vector2d
    { val[0] + [ val[2] ] }

expr:
    expr '+' expr
    { val[0] + val[2] }

  | expr '-' expr
    { val[0] - val[2] }

  | expr '*' expr
    { val[0] * val[2] }

  | expr '/' expr
    { val[0] / val[2] }

  | expr TK_POW expr
    { val[0] ** val[2] }

  | '(' expr ')'
    { val[1] }

  | '-' expr =NEG
    { -val[1] }

  | TK_NUMBER
    { Float(val[0]) }

  | TK_WORD
    { @vars[val[0]] }

  | var_assign

obj:
    TK_BOX point ',' vector
    { Box.new(*val.values_at(1, 3)) }

  | TK_BOX point ',' point
    {
      val[3][0] -= val[1][0]
      val[3][1] -= val[1][1]
      val[3][2] -= val[1][2]
      Box.new(*val.values_at(1, 3))
    }

  | TK_RECT TK_SYMBOL ',' point ',' vector2d
    { Rect.new(*val.values_at(1, 3, 5)) }

  | TK_RECT TK_SYMBOL ',' point ',' point2d
    {
      vec2d = point2d_sub_point(val[1], val[5], val[3])
      Rect.new(val[1], val[3], vec2d)
    }

  | TK_TRIANGLE TK_SYMBOL ',' point ',' vector2d ',' vector2d
    {
      Triangle.new(*val.values_at(1, 3, 5, 7))
    }

  | TK_TRIANGLE TK_SYMBOL ',' point ',' point2d ',' point2d
    {
      vec2d1 = point2d_sub_point(val[1], val[5], val[3])
      vec2d2 = point2d_sub_point(val[1], val[7], val[3])
      Triangle.new(val[1], val[3], vec2d1, vec2d2)
    }

  | TK_CIRCLE TK_SYMBOL ',' point ',' expr
    { Circle.new(*val.values_at(1, 3, 5)) }

  | TK_ELLIPSE TK_SYMBOL ',' point ',' expr ',' expr
    { Ellipse.new(*val.values_at(1, 3, 5, 7)) }

  | TK_POLYGON TK_SYMBOL ',' point ',' vector2d_ary
    { Polygon.new(*val.values_at(1, 3, 5)) }

  | TK_POLYGON TK_SYMBOL ',' point ',' point2d_ary
    {
      ary = []
      val[5].each do |pnt|
        vec2d = point2d_sub_point(val[1], pnt, val[3])
        ary << vec2d
      end
      Polygon.new(val[1], val[3], ary)
    }

  | TK_LINE TK_SYMBOL ',' point ',' vector2d_ary
    { Line.new(*val.values_at(1, 3, 5)) }

  | TK_LINE TK_SYMBOL ',' point ',' point2d_ary
    {
      ary = []
      val[5].each do |pnt|
        vec2d = point2d_sub_point(val[1], pnt, val[3])
        ary << vec2d
      end
      Line.new(val[1], val[3], ary)
    }

  | TK_SWEEP TK_SYMBOL ',' expr ',' obj
    { Sweep.new(*val.values_at(1, 3, 5)) }

  | '[' objs ']'
    {
      val[1]
    }

objs:
    obj
    { Objs.new(val[0]) }

  | objs obj
    {
      val[0] << val[1]
      val[0]
    }


---- header

require 'obj3d/viewer'

class Obj
  attr_accessor :type

  def type=(type)
    @type = type
  end

  def point_add_val(axis, point, val)
    new_point = []
    case axis
    when :X
      new_point = [point[0] + val, point[1]      , point[2]      ]
    when :Y
      new_point = [point[0]      , point[1] + val, point[2]      ]
    when :Z
      new_point = [point[0]      , point[1]      , point[2] + val]
    else
      raise "unknown axis: #{axis}"
    end
    new_point
  end
  private :point_add_val
end

class Rect < Obj
  def initialize(axis, point, vector2d, type = nil)
    @axis = axis
    @point = point
    @vector2d = vector2d
    @type = type
  end

  def draw(viewer, rgb)
    point = @point.map{|v| v * viewer.scale}
    vector2d = @vector2d.map{|v| v * viewer.scale}
    viewer.rect(point, @axis, vector2d, rgb, true)
  end

  def move(axis, val)
    point = point_add_val(axis, @point, val)
    self.class.new(@axis, point, @vector2d, @type)
  end
end

class Triangle < Obj
  def initialize(axis, point, vector2d_a, vector2d_b)
    @polygon = Polygon.new(axis, point, [vector2d_a, vector2d_b])
  end

  def draw(viewer, rgb)
    @polygon.draw(viewer, rgb)
  end

  def move(axis, val)
    polygon = @polygon.move(axis, val)
    polygon.type = @type
    polygon
  end
end

class Ellipse < Obj
  def initialize(axis, point, ru, rv)
    warn 'not yet Ellipse'
  end
end

class Circle < Obj
  def initialize(axis, point, r, type = nil)
    @axis = axis
    @point = point
    @r = r
    @type = type
  end

  def draw(viewer, rgb)
    point = @point.map{|v| v * viewer.scale}
    r = @r * viewer.scale
    viewer.circle(point, @axis, r, rgb, true)
  end

  def move(axis, val)
    point = point_add_val(axis, @point, val)
    self.class.new(@axis, point, @r, @type)
  end
end

class Box < Obj
  def initialize(point, vector)
    @point = point
    @vector = vector
  end

  def draw(viewer, rgb)
    point = @point.map{|v| v * viewer.scale}
    vector = @vector.map{|v| v * viewer.scale}
    viewer.cube(point, vector, rgb, true)
  end
end

class Sweep < Obj
  def initialize(axis, val, obj)
    @obj = obj
    @obj2 = obj.move(axis, val)
  end

  def type=(type)
    super
    @obj.type = @obj2.type = type
  end

  def draw(viewer, rgb)
    @obj.draw(viewer, rgb)
    @obj2.draw(viewer, rgb)
  end
end

class Polygon < Obj
  def initialize(axis, point, vector2d_ary, type = nil)
    @axis = axis
    @point = point
    @vector2d_ary = vector2d_ary
    @type = type
  end

  def draw(viewer, rgb)
    points = [@point]
    @vector2d_ary.each do |v2d|
      case @axis
      when :X
        points << [@point[0]         , @point[1] + v2d[0], @point[2] + v2d[1]]
      when :Y
        points << [@point[0] + v2d[0], @point[1]         , @point[2] + v2d[1]]
      when :Z
        points << [@point[0] + v2d[0], @point[1] + v2d[1], @point[2]         ]
      else
        raise 'bug'
      end
    end
    points = points.map{|point| point.map{|v| v * viewer.scale}}
    viewer.polygon(points, rgb, true)
  end

  def move(axis, val)
    if axis != @axis
      raise 'not implemented'
    end
    point = point_add_val(axis, @point, val)
    self.class.new(@axis, point, @vector2d_ary, @type)
  end
end

class Line < Obj
  def initialize(axis, point, vector2d_ary, type = nil)
    @axis = axis
    @point = point
    @vector2d_ary = vector2d_ary
    @type = type
  end

  def draw(viewer, rgb)
    points = [@point]
    @vector2d_ary.each do |v2d|
      case @axis
      when :X
        points << [@point[0]         , @point[1] + v2d[0], @point[2] + v2d[1]]
      when :Y
        points << [@point[0] + v2d[0], @point[1]         , @point[2] + v2d[1]]
      when :Z
        points << [@point[0] + v2d[0], @point[1] + v2d[1], @point[2]         ]
      else
        raise 'bug'
      end
    end
    points = points.map{|point| point.map{|v| v * viewer.scale}}
    viewer.lines(points, rgb)
  end

  def move(axis, val)
    if axis != @axis
      raise 'not implemented'
    end
    point = point_add_val(axis, @point, val)
    self.class.new(@axis, point, @vector2d_ary, @type)
  end
end

class Objs < Obj
  def initialize(obj)
    @objs = [obj]
  end

  def <<(obj)
    @objs << obj
  end

  def each
    @objs.each do |obj|
      yield obj
    end
  end

  def draw(viewer, rgb)
    @objs.each do |obj|
      obj.draw(viewer, rgb)
    end
  end

  def move(axis, val)
    new_objs = Objs.new(@objs[0].move(axis, val))
    (1 ... @objs.size).each do |i|
      new_objs << @objs[i].move(axis, val)
    end
    new_objs
  end
end


---- inner

NumberPat = '(?: \d*\.\d+(?:[eE][-+]?\d+)? | ' +
  '\d+\.?(?:[eE][-+]?\d+)? )'

attr_reader :obj_ary, :viewer, :min, :max

def initialize
  @viewer = Obj3D::Viewer.new('viewer', 640, 480)
  @min = @max = nil
  @line_no = 1
end

def parse(str)
  @yydebug = false

  @type = nil
  @vars = {}

  @obj_ary = []

  str = str.strip
  @q = []
  until str.empty?
    case str
    when /\A\n/
      @q.push [:TK_EOL, $&]
      str = $'
    when /\A\s+/
      str = $'
    when /\A#.*/
      str = $'
    when /\A---+/
      @q.push [:TK_DASH, $&]
      str = $'
    when /\Aactive\b/
      @q.push [:TK_ACTIVE, $&]
      str = $'
    when /\Abox\b/
      @q.push [:TK_BOX, $&]
      str = $'
    when /\Acircle\b/
      @q.push [:TK_CIRCLE, $&]
      str = $'
    when /\Aellipse\b/
      @q.push [:TK_ELLIPSE, $&]
      str = $'
    when /\Afix\b/
      @q.push [:TK_FIX, $&]
      str = $'
    when /\Aheat\b/
      @q.push [:TK_HEAT, $&]
      str = $'
    when /\Alambda\b/
      @q.push [:TK_LAMBDA, $&]
      str = $'
    when /\Aline\b/
      @q.push [:TK_LINE, $&]
      str = $'
    when /\Ainactive\b/
      @q.push [:TK_INACTIVE, $&]
      str = $'
    when /\Apolygon\b/
      @q.push [:TK_POLYGON, $&]
      str = $'
    when /\Arect\b/
      @q.push [:TK_RECT, $&]
      str = $'
    when /\Asweep\b/
      @q.push [:TK_SWEEP, $&]
      str = $'
    when /\Atriangle\b/
      @q.push [:TK_TRIANGLE, $&]
      str = $'
    when /\Aworld\b/
      @q.push [:TK_WORLD, $&]
      str = $'
    when /\A#{NumberPat}/ox
      @q.push [:TK_NUMBER, $&]
      str = $'
    when /\A[a-zA-Z]\w*/
      @q.push [:TK_WORD, $&]
      str = $'
    when /\A:[a-zA-Z]\w*/
      @q.push [:TK_SYMBOL, eval($&)]
      str = $'
    when /\A\*\*/
      @q.push [:TK_POW, $&]
      str = $'
    else
      c = str[0, 1]
      @q.push [c, c]
      str = str[1 .. -1]
    end
  end
  @q.push [false, '$']   # is optional from Racc 1.3.7
  do_parse
end

def next_token
  token = nil
  loop do
    token = @q.shift
    unless token[0] == :TK_EOL
      break
    end
    @line_no += 1
  end
  token
end

def on_error(error_token_id, error_value, value_stack)
  STDERR.puts "parse error at #{@line_no} on '#{error_value}'"
  exit 1
end

def point2d_sub_point(axis, point2d, point)
  pnt2d = point2d.dup
  case axis
  when :X
    pnt2d[0] -= point[1]
    pnt2d[1] -= point[2]
  when :Y
    pnt2d[0] -= point[0]
    pnt2d[1] -= point[2]
  when :Z
    pnt2d[0] -= point[0]
    pnt2d[1] -= point[1]
  else
    raise "unknown axis '#{axis}'"
  end
  pnt2d
end

def set_min_max(ary)
  ary.each do |v|
    @min = v if @min.nil? or v < @min
    @max = v if @max.nil? or v > @max
  end
end


---- footer

class Obj3D::Viewer
  attr_accessor :scale
end

def draw(viewer, obj_ary, obj_type)
  viewer.delete_list
  obj_ary.each do |obj|
    next unless obj_type == :ALL or obj_type == obj.type
    case obj.type
    when :WORLD
      obj.draw(viewer, [0.5, 0.5, 0.5])
    when :ACTIVE
      obj.draw(viewer, [1, 1, 1])
    when :INACTIVE
      obj.draw(viewer, [0.3, 0.3, 0.3])
    when :FIX
      obj.draw(viewer, [0, 0, 1])
    when :FIXHEAT
      obj.draw(viewer, [1, 0, 0])
    when :HEAT
      obj.draw(viewer, [1, 0, 0])
    when :LAMBDA
      obj.draw(viewer, [0.0, 0.5, 0.0])
    else
      raise "unknown type #{obj.type}"
    end
  end
end

def next_type(type)
  type_ary = [:ALL, :WORLD, :ACTIVE, :INACTIVE, :FIX, :FIXHEAT, :HEAT, :LAMBDA]
  index = type_ary.index(type)
  if index.nil?
    raise 'bug?'
  end
  index = (index + 1) % type_ary.size
  type_ary[index]
end

parser = TcParser.new
parser.parse(ARGF.read)

parser.viewer.scale = 3.0 / (parser.max - parser.min)

obj_type = :ALL

keyboard_proc = Proc.new { |key, x, y|
  case key.chr
  when 'q', 'Q'
    exit
  when 't'
    obj_type = next_type(obj_type)
  when 'w'
    if obj_type != :WORLD
      obj_type = :WORLD
    else
      obj_type = :ALL
    end
  when 'a'
    if obj_type != :ACTIVE
      obj_type = :ACTIVE
    else
      obj_type = :ALL
    end
  when 'A'
    if obj_type != :INACTIVE
      obj_type = :INACTIVE
    else
      obj_type = :ALL
    end
  when 'f'
    if obj_type != :FIX
      obj_type = :FIX
    else
      obj_type = :ALL
    end
  when 'F'
    if obj_type != :FIXHEAT
      obj_type = :FIXHEAT
    else
      obj_type = :ALL
    end
  when 'h'
    if obj_type != :HEAT
      obj_type = :HEAT
    else
      obj_type = :ALL
    end
  when 'l'
    if obj_type != :LAMBDA
      obj_type = :LAMBDA
    else
      obj_type = :ALL
    end
  else
    STDERR.puts "'#{key.chr}' is not assigned"
  end
  draw(parser.viewer, parser.obj_ary, obj_type)
  parser.viewer.newlist
  parser.viewer.redisplay
  STDERR.puts obj_type
}
parser.viewer.keyboard(keyboard_proc)
parser.viewer.light = false

draw(parser.viewer, parser.obj_ary, obj_type)

parser.viewer.mainloop
