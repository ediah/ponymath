use "lib:matrix"
use @mtxmul[Pointer[F64] ref](a: Pointer[F64] tag, b: Pointer[F64] tag,
                  m: USize tag, n: USize tag, k: USize tag)
use @mtxmulpar[None](a: Pointer[F64] tag, b: Pointer[F64] tag,
    c: Pointer[F64] tag, m: USize tag, n: USize tag, k: USize tag,
    i: USize tag, threads: USize tag, semaphore: Pointer[F64] tag)

use "collections"
use "random"

class Matrix
  let m: USize
  let n: USize
  var data: Array[F64]

  new create(m': USize, n': USize) =>
    m = m'
    n = n'
    data = Array[F64](m * n)

  fun ref init(seed: U64, max: F64): Matrix =>
    var random = Rand(seed)
    random.real() // skip first number
    for i in Range(0, m * n) do
      data.push(random.real() * max)
    end
    this
/*
  fun string(): String =>
    var str = "Matrix[" + m.string() + " x " + n.string() + "]\n"

    for i in Range(0, m) do
      for j in Range(0, n - 1) do
        str.append(this(i, j).string() + ", ")
      end
      str.append(this(i, n).string() + "\n")
    end
    str
*/

  fun string(): String =>
    var str = "Matrix[" + m.string() + "x" + n.string() + "]\n"
    let dv = data.values()
    str.append("{")
    for i in Range(0, m - 1) do
      str.append("{")
      for j in Range(0, n - 1) do
        try str.append(dv.next()?.string() + ", ") end
      end
      try str.append(dv.next()?.string() + "},\n") end
    end
    str.append("{")
    for j in Range(0, n - 1) do
      try str.append(dv.next()?.string() + ", ") end
    end
    try str.append(dv.next()?.string() + "}") end
    str.append("}")
    str

  fun apply(i: USize, j: USize): F64 =>
    try data((i * n) + j)? else F64(0) end
    //F64.from[USize](data.size())

  fun mul_naive(y: Matrix): Matrix =>
    var result = Matrix(m, y.n)
    var sum: F64 = 0
    for i in Range(0, m) do
      for j in Range(0, y.n) do
        for k in Range(0, n) do
          sum = sum + (this(i, k) * y(k, j))
        end
        result.data.push(sum)
        sum = 0
      end
    end
    result

  fun transpose(): Matrix =>
    var res = Matrix(m, n)
    for i in Range(0, m) do
      for j in Range(0, n) do
        res.data.push(this(j, i))
      end
    end
    res

  fun mul_transposed(yT: Matrix): Matrix =>
    var c = Matrix(m, yT.m)
    var cpoint = @mtxmul(this.data.cpointer(),
                           yT.data.cpointer(),
                         this.m, this.n, yT.m)
    c.data = Array[F64].from_cpointer(cpoint, m * yT.m)
    c

  fun ref initzero(): Matrix =>
    for i in Range(0, m * n) do
      data.push(0)
    end
    this

  fun iszero(): Bool =>
    for i in Range(0, m) do
      for j in Range(0, n) do
        if this(i, j) == 0 then
          return true
        end
      end
    end
    false

  fun mul_actors(yT: Matrix): Matrix =>
    let threads: USize = 2
    var c = Matrix(m, yT.m).initzero()
    var sem = Matrix(1, threads - 1).initzero()
    var marray = Array[MatrixActor](threads - 1)
    for i in Range(1, threads) do
      marray.push(MatrixActor(
        this.data.cpointer(), // Указатель на начало A
        yT.data.cpointer(),   // Указатель на начало B
        c.data.cpointer(),    // Указатель на начало C
        this.m,               // Кол-во строк A
        this.n,               // Кол-во столбцов A = кол-во строк B
        yT.m,                 // Кол-во столбцов B
        i,                    // Номер части
        threads,              // Сколько всего частей
        sem.data.cpointer()   // Семафор
        ))
    end

    for i in Range(1, threads + 1) do
      try marray(i - 1)?.multiplier() end
    end

    @mtxmulpar(
      this.data.cpointer(), // Указатель на начало A
      yT.data.cpointer(),   // Указатель на начало B
      c.data.cpointer(),    // Указатель на начало C
      this.m,               // Кол-во строк A
      this.n,               // Кол-во столбцов A = кол-во строк B
      yT.m,                 // Кол-во столбцов B
      threads,              // Номер части
      threads,              // Сколько всего частей
      sem.data.cpointer()   // Семафор
      )

    while sem.iszero() do
      dummy()
    end

    c


  fun dummy() =>
    None

  fun mul(y: Matrix): Matrix =>
    this.mul_actors(y.transpose())

actor MatrixActor
  var a: Pointer[F64] tag
  var b: Pointer[F64] tag
  var c: Pointer[F64] tag
  var m: USize tag
  var n: USize tag
  var k: USize tag
  var i: USize tag
  var threads: USize tag
  var semaphore: Pointer[F64] tag

  new create(a': Pointer[F64] tag, b': Pointer[F64] tag, c': Pointer[F64] tag,
             m': USize tag, n': USize tag, k': USize tag, i': USize tag,
             threads': USize tag, semaphore': Pointer[F64] tag) =>
    a = a'
    b = b'
    c = c'
    m = m'
    n = n'
    k = k'
    i = i'
    threads = threads'
    semaphore = semaphore'

  be multiplier() =>
    @mtxmulpar(a, b, c, m, n, k, i, threads, semaphore)
