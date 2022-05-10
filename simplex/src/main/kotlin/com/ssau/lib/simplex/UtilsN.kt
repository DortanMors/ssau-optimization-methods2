package com.ssau.lib.simplex

import org.jetbrains.kotlinx.multik.api.d2arrayIndices
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray
import org.jetbrains.kotlinx.multik.api.zeros
import org.jetbrains.kotlinx.multik.ndarray.data.D1Array
import org.jetbrains.kotlinx.multik.ndarray.data.D2Array
import org.jetbrains.kotlinx.multik.ndarray.data.get
import org.jetbrains.kotlinx.multik.ndarray.operations.*

internal typealias Vector = D1Array<Double>
internal typealias Matrix = D2Array<Double>

fun vector(size: Int) = mk.zeros<Double>(size)

fun vector(vararg numbers: Double) = mk.ndarray(numbers)

fun vector(v: List<Double>) = mk.ndarray(v)

fun matrix(vararg vectors: List<Double>): Matrix = mk.ndarray(vectors.toList())

fun matrix(a: List<List<Double>>) = mk.ndarray(a)
fun Vector.asMatrix() = mk.d2arrayIndices(size, size) { row, col ->
    if (row == 0) this[col]
    else 0.0
}

fun Matrix.asVector() = transpose()[0]

fun Matrix.rows() = shape[0]

fun Matrix.cols() = shape[1]

operator fun Matrix.times(right: Vector): Matrix = this * right.asMatrix()

operator fun Vector.minus(right: Matrix): Vector = this - right.asVector()

fun Matrix.withRow(vector: Vector) = append(vector).reshape(shape[0]+1, shape[1])

fun Matrix.withColumn(vector: Vector) = withColumn(vector.toList())

fun Matrix.withColumn(list: List<Double>) = mk.d2arrayIndices(shape[0], shape[1] + 1) { i, j ->
    if (j < shape[1]) this[i][j] else list[i]
}
