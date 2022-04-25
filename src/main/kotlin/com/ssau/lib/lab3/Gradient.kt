package com.ssau.lib.lab3

import com.ssau.lib.*
import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import org.jetbrains.kotlinx.multik.ndarray.operations.plus
import org.jetbrains.kotlinx.multik.ndarray.operations.times


object Gradient : MinFinderN {
    override fun findMinN(f: (x: Vector) -> Double, x0: Vector, x1: Vector, eps: Double, maxIterations: Int): Vector {
        var xI = x0.copy()
        var xI1 = x0.copy()
        var cntr = 0
        while (true) {
            cntr++
            if (cntr == maxIterations) {
                break
            }
            xI1 = xI - gradient(f, xI, eps)
            xI1 = DichotomyN.findMinN(f, xI, xI1, eps, maxIterations)
            if ((xI1 - xI).abs() < eps) break
            xI = xI1
        }
        Logger.logIters(javaClass.simpleName, cntr)
        return (xI1 + xI) * 0.5
    }
}