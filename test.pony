actor Main
  new create(env: Env) =>
    var n: USize = 2000
    var mtx1 = Matrix(n, n).init(200, 1)
    var mtx2 = Matrix(n, n).init(4000, 1)

    //env.out.print(mtx1.string())
    //env.out.print(mtx2.string())

    var mtx3 = mtx1 * mtx2

    //env.out.print(mtx3.string())
