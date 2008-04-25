class TCCParser

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
      @model = Model.new(val[2], val[4], val.values_at(6, 8, 10))
    }

  | TK_DASH TK_WORLD point ',' point ',' expr ',' expr ',' expr
    {
      val[4][0] -= val[2][0]
      val[4][1] -= val[2][1]
      val[4][2] -= val[2][2]
      @model = Model.new(val[2], val[4], val.values_at(6, 8, 10))
    }

commands:
    command

  | commands command

command:
    TK_DASH TK_ACTIVE obj
    {
      @model.active(val[2])
    }

  | TK_DASH TK_INACTIVE obj
    {
      @model.inactive(val[2])
    }

  | TK_DASH TK_FIX expr obj
    {
      @model.fix(val[2], val[3])
    }

  | TK_DASH TK_FIX TK_HEAT expr obj
    {
      @model.fixheat(val[3], val[4])
    }

  | TK_DASH TK_HEAT TK_FIX expr obj
    {
      @model.fixheat(val[3], val[4])
    }

  | TK_DASH TK_HEAT expr obj
    {
      @model.heat(val[2], val[3])
    }

  | TK_DASH TK_LAMBDA expr obj
    {
      @model.lambda(val[2], val[3])
    }

point:
    '(' expr ',' expr ',' expr ')'
    {
      point = val.values_at(1, 3, 5)
      point
    }

point2d:
    '(' expr ',' expr ')'
    {
      point2d = val.values_at(1, 3)
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
      vector
    }

vector2d:
    '<' expr ',' expr '>'
    {
      vector2d = val.values_at(1, 3)
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

  | TK_PI

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

  | TK_CIRCLEPERI TK_SYMBOL ',' point ',' expr ',' expr ',' expr
    { CirclePeri.new(*val.values_at(1, 3, 5, 7, 9)) }

  | TK_ELLIPSE TK_SYMBOL ',' point ',' expr ',' expr
    { Ellipse.new(*val.values_at(1, 3, 5, 7)) }

  | TK_ELLIPSEPERI TK_SYMBOL ',' point ',' expr ',' expr ',' expr ',' expr
    { EllipsePeri.new(*val.values_at(1, 3, 5, 7, 9, 11)) }

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
    { Objs[val[0]] }

  | objs obj
    {
      val[0] << val[1]
      val[0]
    }


---- header

require 'tcc'

include TCC


---- inner

NumberPat = '(?: \d*\.\d+(?:[eE][-+]?\d+)? | ' +
  '\d+\.?(?:[eE][-+]?\d+)? )'

def initialize
  @model = nil
  @vars = {}
  @line_no = 1
end

def parse(str)
  @yydebug = false

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
    when /\Acircleperi\b/
      @q.push [:TK_CIRCLEPERI, $&]
      str = $'
    when /\Aellipse\b/
      @q.push [:TK_ELLIPSE, $&]
      str = $'
    when /\Aellipseperi\b/
      @q.push [:TK_ELLIPSEPERI, $&]
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
    when /\Api\b/
      @q.push [:TK_PI, Math::PI]
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


---- footer

parser = TCCParser.new
parser.parse(ARGF.read)
