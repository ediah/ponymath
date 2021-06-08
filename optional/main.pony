actor Main
  new create(env: Env) =>
    var n: USize = 2000
    var a = Matrix(n, n).init(200, 1)
    var b = Matrix(n, n).init(4000, 1)

    var c = a * b
