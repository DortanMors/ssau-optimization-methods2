package com.ssau.lib

interface MinFinderN {
    fun findMinN(f: (x: Vector) -> Double, x0: Vector, x1: Vector, eps: Double, maxIterations: Int): Vector
}