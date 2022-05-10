package com.ssau.lib.simplex

import org.junit.jupiter.api.BeforeEach
import org.junit.jupiter.api.Test

internal class SimplexTest {

    private lateinit var a: List<List<Double>>
    private lateinit var b: List<Double>
    private lateinit var c: List<Double>
    private lateinit var c1: List<Double>
    private lateinit var c2: List<Double>
    private lateinit var c3: List<Double>
    private lateinit var signsLess: List<Sign>
    private lateinit var signsGreater: List<Sign>
    private lateinit var signsEqual: List<Sign>

    @BeforeEach
    fun setUp() {
        a = listOf(
            listOf(-2.0, 6.0),
            listOf(3.0, 2.0),
            listOf(2.0, -1.0)
        )
        b = listOf(40.0, 28.0, 14.0)
        c = listOf(2.0, 3.0)
        c1 = listOf(-2.0, 3.0)
        c2 = listOf(2.0, 1.0)
        c3 = listOf(-2.0, -1.0)
        signsLess = listOf(Sign.Less, Sign.Less, Sign.Less)
        signsGreater = listOf(Sign.Greater, Sign.Greater, Sign.Greater)
        signsEqual = listOf(Sign.Equal, Sign.Equal, Sign.Equal)
    }

    @Test
    fun solve() {
        println(" f(x,c) =  2x1 + 3x2;\n arg_max = {4, 8}, f(arg_max) = 32")
        println(" |-2x1 + 6x2 <= 40")
        println(" | 3x1 + 2x2 <= 28")
        println(" | 2x1 -  x2 <= 14\n")
        val simplex = Simplex(a, c, signsLess, b)
        simplex.findMax().solve()
        println("\n f(x,c) = -2x1 + 3x2;\n arg_min = {7, 0}, f(arg_min) =-14\n")
        simplex.withPrices(c1).findMin().solve()
        println("/////////////////////////////")
        println(" f(x,c) =  2x1 + 3x2;\n arg_min = {62/5, 54/5}, f(arg_max) = 57 1/5")
        println(" |-2x1 + 6x2 >= 40")
        println(" | 3x1 + 2x2 >= 28")
        println(" | 2x1 -  x2 >= 14\n")
        simplex.withPrices(c2).withInequalities(signsGreater).findMin().solve()
        println("/////////////////////////////")
        println(" f(x,c) =  -2x1 - x2;\n arg_min = {62/5, 54/5}, f(arg_max) = -35 3/5")
        simplex.withPrices(c3).findMax().solve()
        println("/////////////////////////////")
        println(" f(x,c) =  2x1 + 3x2;\n arg_min = {none, none}, f(arg_max) = none")
        simplex.withPrices(c).withInequalities(signsEqual).findMax().solve()
    }
}