package com.ssau.lib

import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import org.jetbrains.kotlinx.multik.ndarray.operations.plus
import org.jetbrains.kotlinx.multik.ndarray.operations.times

object DichotomyN : MinFinderN {
    override fun findMinN(f: (x: Vector) -> Double, x0: Vector, x1: Vector, eps: Double, maxIterations: Int)
    : Vector {
        var start = x0
        var end = x1
        val dir = calcDir(start, end)
        var i = 0
        while (i < maxIterations) {
            if ((end - start).abs() < eps) break
            val xc = (start + end) * 0.5
            if (f(xc + dir * eps) > f(xc - dir * eps))
                end = xc
            else
                start = xc
            ++i
        }
        Logger.logIters(javaClass.simpleName, i)
        return (start + end) * 0.5
    }
}