package com.ssau.lib.lab3

import com.ssau.lib.*
import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import org.jetbrains.kotlinx.multik.ndarray.operations.plus
import org.jetbrains.kotlinx.multik.ndarray.operations.times
import kotlin.math.pow


object ConjGradient : MinFinderN {
    override fun findMinN(f: (x: Vector) -> Double, x0: Vector, x1: Vector, eps: Double, maxIterations: Int): Vector {
        var xI: Vector = x0.copy()

        var xI1: Vector = x0.copy()

        val sI: Vector = gradient(f, x0, eps) * (-1.0)
        var sI1: Vector

        var omega: Double

        var cntr = 0

        while (true) {
            cntr++
            if (cntr == maxIterations) {
                break
            }
            xI1 = xI + sI
            xI1 = DichotomyN.findMinN(f, xI, xI1, eps, maxIterations)
            if ((xI1 - xI).abs() < eps) {
                break
            }
            sI1 = gradient(f, xI1, eps)
            omega = sI1.abs().pow(2) / sI.abs().pow(2)
            sI * omega - sI1
            xI = xI1
        }

        Logger.logIters(javaClass.name, cntr)
        return (xI1 + xI) * 0.5
    }
}