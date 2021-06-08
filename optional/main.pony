actor Main
  new create(env: Env) =>
    var n: USize = 2000
    var a = Matrix(n, n).init(200, 1)
    var b = Matrix(n, n).init(4000, 1)

    //env.out.print(a.string())
    //env.out.print(b.string())

    var c = a * b

    //env.out.print(c.string())
    //env.out.print(a.mul_naive(b).string())
