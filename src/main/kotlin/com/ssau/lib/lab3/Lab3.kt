package com.ssau.lib.lab3

import com.ssau.lib.Vector
import com.ssau.lib.abs
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray
import org.jetbrains.kotlinx.multik.ndarray.data.get
import org.jetbrains.kotlinx.multik.ndarray.operations.minus
import kotlin.math.pow

class Lab3 {

    companion object {
        @JvmStatic
        fun main(args: Array<String>) {
            val minX: Vector = mk.ndarray(mk[2.0, -3.0])
            val methods = listOf(Gradient, ConjGradient)
            methods.forEach { method ->
                val resultX = method.findMinN(
                    { v: Vector -> (v[0] - minX[0]).pow(2) + (v[1] - minX[1]).pow(2) },
                    mk.ndarray(mk[-1.5, -1.5]),
                    mk.ndarray(mk[1.0, 1.0]),
                    1e-6,
                    1000
                )
                println("result $resultX error = ${(resultX - minX).abs()} by ${ method.javaClass.simpleName }")
            }
        }
    }
}