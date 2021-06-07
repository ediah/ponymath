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
  let _m: USize
  let _n: USize
  let _size: USize
  var data: Array[F64]
  var _aligned: Bool = false
  var _collected: Bool = false

  new create(m': USize, n': USize) =>
    _m = m'
    _n = n'
    _size = m' * n'
    data = Array[F64](_size)

  fun size(): USize =>
    _size

  fun ref init(seed: U64, max: F64): Matrix =>
    var random = Rand(seed)
    random.real() // skip first number
    for i in Range(0, _size) do
      data.push(random.real() * max)
    end
    this

  fun ref update(i: USize, value: F64) =>
    try data.update(i, value)? end

  fun string(): String =>
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

  fun apply(i: USize, j: USize): F64 =>
    try data((i * _n) + j)? else F64(0) end

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

  fun transpose(): Matrix =>
    var res = Matrix(_n, _m)
    for i in Range(0, _m) do
      for j in Range(0, _n) do
        res.data.push(this(j, i))
      end
    end
    res

  fun copy(): Matrix =>
    var res = Matrix(_n, _m)
    for i in Range(0, _m) do
      for j in Range(0, _n) do
        res.data.push(this(i, j))
      end
    end
    res

  fun mul_transposed(yT: Matrix): Matrix =>
    var c = Matrix(_m, yT._m)
    var cpoint = @mtxmul(
      this.data.cpointer(), // A
      yT.data.cpointer(),   // B
      _m,    // A rows
      _n,    // A columns
      yT._m  // B rows
    )
    c.data = Array[F64].from_cpointer(cpoint, yT._size)
    c

  fun ref initzero(): Matrix =>
    for i in Range(0, _size) do
      data.push(0)
    end
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

  fun ref realign(): Matrix =>
    data.from_cpointer(
      @realign(data.cpointer(), _size),
      _size
    )
    _aligned = true
    this

  fun ref mul_actors(y: Matrix): Matrix =>
    let cores: USize = 2
    var res: Matrix iso^ = recover res.create(_m, y._n).initzero() end

    if not _aligned then
      realign()
    end

    var block: USize = y._n / cores
    var collector = MatrixCollector(res)
    var workers = Array[MatrixActor](cores)

    for i in Range(0, cores) do
      var k: USize = block
      if i == (cores - 1) then
        k = block + (y._n %% cores)
      end
      var from: USize = block * i
      var worker = MatrixActor(
        copy().data.cpointer(),
        y.sliceT(from, from + k)
         .realign()
         .data.cpointer(),
        _m, _n, k,
        collector
      )
      collector.add(worker, from * _m)
      workers.push(worker)
    end

    for i in Range(0, cores) do
      try workers(i)?.multiplier() end
    end

    _wait_for_collection(collector, res)

    res


  fun ref collection_done() =>
    _collected = true

  fun _wait_for_collection(collector: MatrixCollector, target: Matrix iso^) =>
    while not target._collected do
      collector.answer(target)
    end

  fun ref mul(y: Matrix): Matrix =>
    this.mul_actors(y)

  new from_array(arr: Array[F64], m: USize, n: USize) =>
    data = arr
    _m = m
    _n = n
    _size = m * n

actor MatrixCollector
  var _workers: Array[MatrixActor] = _workers.create()
  var _shifts: Array[USize] = _shifts.create()
  let _result: Matrix iso

  new create(result: Matrix iso^) =>
    _result = result

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

      _workers.delete(nworker)?
      _shifts.delete(nworker)?
    end

  be answer(done: Matrix iso) =>
    if _workers.size() == 0 then
      done.collection_done()
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
