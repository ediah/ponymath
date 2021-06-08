use "lib:matrix"

// Matrix multiplication
use @mtxmul[Pointer[F64] ref](
  a: Pointer[F64] tag, // Left matrix
  b: Pointer[F64] tag, // Right matrix (transposed)
  m: USize tag,        // Rows in left matrix
  n: USize tag,        // Columns in left matrix
  k: USize tag         // Rows in transposed right matrix
)

// Make matrix data be aligned
use @realign[Pointer[F64] ref](p: Pointer[F64] tag, size: USize tag)
use @debug_int[None](i: USize)
use @debug_point[None](p: Pointer[F64] tag)

use "collections"
use "random"
use "time"

class Matrix
  let _m: USize                 // Rows num
  let _n: USize                 // Columns num
  let _size: USize              // Elements num
  var data: Array[F64]          // Data
  var _collected: Bool = false  // Is it collected from actors
  let _eps: F64 = 0.000001      // Trust interval for approximations

  new create(m': USize, n': USize) =>
    _m = m'
    _n = n'
    _size = m' * n'
    data = Array[F64](_size)

  new from_array(arr: Array[F64], m: USize, n: USize) =>
    data = arr
    _m = m
    _n = n
    _size = m * n

  /* Standard functions
   */
  fun size(): USize =>
    _size

  fun apply(i: USize, j: USize): F64 =>
    try data((i * _n) + j)? else F64(0) end

  fun ref update(i: USize, value: F64) =>
    try data.update(i, value)? end

  fun string(): String iso^ =>
    """
    This function gives matrix representation which
    is suitable for copy-pasting it in Wolfram Alpha
    """
    var str = "Matrix[" + _m.string() + "x" + _n.string() + "]\n"
    let dv = data.values()
    str.append("{")
    for i in Range(0, _m - 1) do
      str.append("{")
      for j in Range(0, _n - 1) do
        try str.append(dv.next()?.string() + ", ") end
      end
      try str.append(dv.next()?.string() + "},\n") end
    end
    str.append("{")
    for j in Range(0, _n - 1) do
      try str.append(dv.next()?.string() + ", ") end
    end
    try str.append(dv.next()?.string() + "}") end
    str.append("}")
    str

  fun ref init(seed: U64, max: F64): Matrix =>
    var random = Rand(seed)
    random.real() // skip first number
    for i in Range(0, _size) do
      data.push(random.real() * max)
    end
    this

  fun ref initzero(): Matrix =>
    for i in Range(0, _size) do
      data.push(0)
    end
    this

  fun copy(): Matrix =>
    var res = Matrix(_n, _m)
    for i in Range(0, _m) do
      for j in Range(0, _n) do
        res.data.push(this(i, j))
      end
    end
    res

  fun transpose(): Matrix =>
    var res = Matrix(_n, _m)
    for i in Range(0, _n) do
      for j in Range(0, _m) do
        res.data.push(this(j, i))
      end
    end
    res

  fun ref realign(): Matrix =>
    data.from_cpointer(
      @realign(data.cpointer(), _size),
      _size
    )
    this

  fun slice(from: USize, to: USize): Matrix =>
    """
    Make slice of matrix by rows including "to"
    """
    // Check for negative size?
    var res = Matrix(to - from, _n)
    for i in Range(from, to + 1) do
      for j in Range(0, _n) do
        res.data.push(this(i, j))
      end
    end
    res

  fun sliceT(from: USize, to: USize): Matrix =>
    """
    Make slice of matrix by columns including "to"
    and transpose it
    """
    // Check for negative size?
    var res = Matrix(to - from, _m)
    for i in Range(from, to + 1) do
      for j in Range(0, _m) do
        res.data.push(this(j, i))
      end
    end
    res

  /* Multiplication functions
   * mul_naive      -- pure Pony                                    ~0.5 GFLOPS
   * mul_transposed -- C ABI, transposed, AVX/AVX2, cache-optimized   ~5 GFLOPS
   * mul_actors     -- same + multithreaded            (on 2 cores)  ~10 GFLOPS
   *
   * By default, mul is mul_actors.
   */
  fun mul_naive(y: Matrix): Matrix =>
    var result = Matrix(_m, y._n)
    var sum: F64 = 0
    for i in Range(0, _m) do
      for j in Range(0, y._n) do
        for k in Range(0, _n) do
          sum = sum + (this(i, k) * y(k, j))
        end
        result.data.push(sum)
        sum = 0
      end
    end
    result

  fun mul_transposed(yT: Matrix): Matrix =>
    var c = Matrix(_m, yT._m)
    var cpoint = @mtxmul(
      this.data.cpointer(),
      yT.data.cpointer(),
      _m, _n, yT._m
    )
    c.data = Array[F64].from_cpointer(cpoint, yT._size)
    c

  fun ref mul_actors(y: Matrix): Matrix =>
    let cores: USize = 2
    var res: Matrix iso^ = recover res.create(_m, y._n).initzero() end
    var yT = y.transpose()

    var block: USize = _m / cores
    var collector = MatrixCollector(res, cores)
    var workers = Array[MatrixActor](cores)
    var k: USize = block

    for i in Range(0, cores - 1) do
      var from: USize = block * i
      var worker = MatrixActor(
        slice(from, from + k)
         .data.cpointer(),
        yT.data.cpointer(),
        k, _n, yT._m,
        collector
      )
      collector.add(worker, from * yT._m)
      workers.push(worker)
    end

    for i in Range(0, cores - 1) do
      try workers(i)?.multiplier() end
    end

    k = block + (_m %% cores)
    var point: Pointer[F64] = @mtxmul(
      slice(block * (cores - 1), _m)
      .data.cpointer(),
      yT.data.cpointer(),
      k, _n, yT._m
    )

    try
      var arr = Array[F64].from_cpointer(point, k * yT._m)
      var shift = (res._n * block) * (cores - 1)
      for i in Range(0, k * yT._m) do
        res(shift + i) = arr(i)?
      end
    end

    _wait_for_collection(collector, res)

    res

  fun ref mul(y: Matrix): Matrix =>
    this.mul_transposed(y)

  /* Logic functions
   */
  fun box eq(y: Matrix box): Bool val =>
    if (this._m != y._m) or (this._n != y._n) then
      return false
    end
    for i in Range(0, _m) do
      for j in Range(0, _n) do
        if (this(i, j) - y(i, j)) > _eps then
          return false
        end
      end
    end
    true

  fun box ne(y: Matrix box): Bool val =>
    not (this == y)



  /* Misc functions. They are for parallelism.
   */
  fun ref collection_done() =>
    _collected = true

  fun _wait_for_collection(collector: MatrixCollector, target: Matrix iso^) =>
    while not target._collected do
      collector.answer({()(done = target) => done.collection_done()})
    end

actor MatrixCollector
  var _workers: Array[MatrixActor] = _workers.create()
  var _shifts: Array[USize] = _shifts.create()
  let _result: Matrix iso
  let _cores: USize

  new create(result: Matrix iso^, cores: USize) =>
    _result = result
    _cores = cores

  be add(worker: MatrixActor, shift: USize) =>
    _workers.push(worker)
    _shifts.push(shift)

  be collect(worker: MatrixActor, p: Pointer[F64] iso^, m: USize, n: USize) =>
    var arr: Array[F64] = arr.from_cpointer(p, m * n)

    try
      var nworker = _workers.find(worker)?
      for i in Range(0, m * n) do
        _result(_shifts(nworker)? + i) = arr(i)?
      end

      _shifts.delete(nworker)?
      _workers.delete(nworker)?
    end

  be answer(callback: {ref (): None} iso) =>
    if _workers.size() == 0 then
      callback()
    end

actor MatrixActor
  var _a: Pointer[F64] tag // A matrix, left
  var _b: Pointer[F64] tag // B matrix, right
  var _m: USize val // A rows
  var _n: USize val // A columns
  var _k: USize val // B columns if it's not transposed
  let _collector: MatrixCollector

  new create(a': Pointer[F64] tag, b': Pointer[F64] tag,
             m': USize val, n': USize val, k': USize val,
             collector: MatrixCollector) =>
    _a = a'
    _b = b'
    _m = m'
    _n = n'
    _k = k'
    _collector = collector

  be multiplier() =>
    var point: Pointer[F64] iso^ = recover @mtxmul(_a, _b, _m, _n, _k) end
    _collector.collect(this, point, _m, _k)
