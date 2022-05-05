package com.ssau.lib.lab4

import com.ssau.lib.Vector
import com.ssau.lib.abs
import org.jetbrains.kotlinx.multik.api.mk
import org.jetbrains.kotlinx.multik.api.ndarray
import org.jetbrains.kotlinx.multik.ndarray.data.get
import org.jetbrains.kotlinx.multik.ndarray.operations.minus

class Lab4 {
    companion object {
        @JvmStatic
        fun main(args: Array<String>) {
            val minX: Vector = mk.ndarray(mk[2.0, 2.0])
            val methods = listOf(NewtoneRaphson)

            methods.forEach { method ->
                val resultX = method.findMinN(
                        { x: Vector -> (x[0] - minX[0]) * x[0] + (x[1] - minX[1]) * x[1]},
                    mk.ndarray(mk[0.0, 0.0]),
                    0.001,
                    1000
                )
                println("result $resultX error = ${(resultX - minX).abs()} by ${ method.javaClass.simpleName }")
            }
        }
    }
}