package com.ssau.lib

import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.zeros
import org.jetbrains.kotlinx.multik.ndarray.data.D1Array
import org.jetbrains.kotlinx.multik.ndarray.data.get
import org.jetbrains.kotlinx.multik.ndarray.data.set
import org.jetbrains.kotlinx.multik.ndarray.operations.map
import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import org.jetbrains.kotlinx.multik.ndarray.operations.sum
import kotlin.math.sqrt

typealias Vector = D1Array<Double>

fun Vector.abs() = sqrt(map { it * it }.sum())

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

fun calcDir(a: Vector, b: Vector) = (b - a).normalized()

fun Vector.normalized() = (1 / abs()).let { norm -> map { it * norm } }