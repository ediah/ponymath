use "ponytest"

actor Main is TestList
  new create(env: Env) =>
    PonyTest(env, this)

  new make() =>
    None

  fun tag tests(test: PonyTest) =>
    test(_TestMul)

class iso _TestMul is UnitTest
  fun name(): String => "multiplication"

  fun apply(h: TestHelper) =>
    var a = Matrix(3, 7).init(200, 1)
    var b = Matrix(7, 10).init(4000, 1)
    h.assert_eq[Matrix](a.mul_naive(b), a * b)
    a = Matrix(100, 8).init(200, 1)
    b = Matrix(8, 100).init(4000, 1)
    h.assert_eq[Matrix](a.mul_naive(b), a * b)
    a = Matrix(500, 500).init(200, 1)
    b = Matrix(500, 500).init(4000, 1)
    h.assert_eq[Matrix](a.mul_naive(b), a * b)
