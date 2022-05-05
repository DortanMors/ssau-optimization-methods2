package com.ssau.lib.lab4

import com.ssau.lib.*
import org.jetbrains.kotlinx.multik.api.linalg.dot
import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import org.jetbrains.kotlinx.multik.ndarray.operations.plus
import org.jetbrains.kotlinx.multik.ndarray.operations.times


object NewtoneRaphson : MinFinderN {
    override fun findMinN(f: (Vector) -> Double, x0: Vector, x1: Vector, eps: Double, maxIterations: Int): Vector {
        var xI: Vector = Vector(x0)
        var xI1: Vector = Vector(x0)
        var cntr = -1
        while (true) {
            cntr++
            if (cntr == maxIterations) {
                break
            }
            //val hess = hessian(f, xI, eps);
            //val grad = gradient(f, xI, eps);
            xI1 = xI - ((hessian(f, xI, eps)).invert() dot gradient(f, xI, eps))

            if ((xI1 - xI).abs() < eps)
                break
            xI = xI1
        }
        Logger.logIters(javaClass.simpleName, cntr)
        return (xI1 + xI) * 0.5
    }

    fun findMinN(f: (Vector) -> Double, x0: Vector, eps: Double, maxIterations: Int): Vector =
        findMinN(f, x0, x0, eps, maxIterations)
}
