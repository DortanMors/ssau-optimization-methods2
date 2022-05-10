package com.ssau.lib

import org.jetbrains.kotlinx.multik.api.d2arrayIndices
import org.jetbrains.kotlinx.multik.api.linalg.inv
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.zeros
import org.jetbrains.kotlinx.multik.ndarray.data.D1Array
import org.jetbrains.kotlinx.multik.ndarray.data.D2Array
import org.jetbrains.kotlinx.multik.ndarray.data.get
import org.jetbrains.kotlinx.multik.ndarray.data.set
import org.jetbrains.kotlinx.multik.ndarray.operations.map
import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import org.jetbrains.kotlinx.multik.ndarray.operations.sum
import org.jetbrains.kotlinx.multik.ndarray.operations.times
import kotlin.math.sqrt


typealias Vector = D1Array<Double>
typealias Matrix = D2Array<Double>

fun Vector.abs() = sqrt(map { it * it }.sum())

fun Vector(value: Vector) = value.copy()

fun Vector(size: Int) = mk.zeros<Double>(size)
fun gradient(f: (Vector) -> Double, x: Vector, eps: Double): Vector {
    val xL: Vector = x.copy()
    val xR: Vector = x.copy()
    val df: Vector = mk.zeros(x.size)
    for (i in x.indices) {
        xL[i] -= eps
        xR[i] += eps
        df[i] = (f.invoke(xR) - f.invoke(xL)) * (0.5 / eps)
        xL[i] += eps
        xR[i] -= eps
    }
    return df
}

fun Vector.asMatrix() = mk.d2arrayIndices(size, size) { row, col ->
    if (row == 0) this[col]
    else 0.0
}

fun Matrix.asVector() = transpose()[0]

fun Matrix.rows() = shape[0]

fun Matrix.cols() = shape[1]

operator fun Matrix.times(right: Vector): Matrix = this * right.asMatrix()

operator fun Vector.minus(right: Matrix): Vector = this - right.asVector()

fun hessian(f: (Vector) -> Double, x: Vector, eps: Double): Matrix =
    mk.d2arrayIndices(x.size, x.size) { row, col ->
        partial2(f, x, row, col, eps)
    }

fun partial2(func: (Vector) -> Double, x: Vector, index_1: Int, index_2: Int, eps: Double): Double {
    if (index_2 !in 0 until x.size) {
        throw IndexOutOfBoundsException("Partial derivative index out of bounds!")
    }
    x[index_2] -= eps
    val fL = partial(func, x, index_1, eps)
    x[index_2] += 2.0 * eps
    val fR = partial(func, x, index_1, eps)
    x[index_2] -= eps
    return (fR - fL) / eps * 0.5
}

fun partial(f: (Vector) -> Double, x: Vector, index: Int, eps: Double): Double {
    if (index !in 0 until x.size) {
        throw IndexOutOfBoundsException("Partial derivative index out of bounds!")
    }
    x[index] += eps
    val fR = f.invoke(x)
    x[index] -= 2.0 * eps
    val fL = f.invoke(x)
    x[index] += eps
    return (fR - fL) / eps * 0.5
}

fun calcDir(a: Vector, b: Vector) = (b - a).normalized()

fun Vector.normalized() = (1 / abs()).let { norm -> map { it * norm } }

fun Matrix.invert(): Matrix = mk.linalg.inv(this)