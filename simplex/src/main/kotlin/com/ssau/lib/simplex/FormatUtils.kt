package com.ssau.lib.simplex

import org.jetbrains.kotlinx.multik.ndarray.data.get
import kotlin.math.abs

fun decimalToRational(value: Double, maxDen: Int = 1000): IntArray {
    var m00: Long = 1
    var m01: Long = 0
    var m10: Long = 0
    var m11: Long = 1
    val number = IntArray(3)
    var ai: Long
    var x: Double
    val sign = if (value >= 0) 1 else -1
    x = abs(value)
    var t: Long
    while (m10 * x.toLong().also { ai = it } + m11 <= maxDen) {
        t = m00 * ai + m01
        m01 = m00
        m00 = t
        t = m10 * ai + m11
        m11 = m10
        m10 = t
        if (x == ai.toDouble()) {
            break
        }
        x = 1 / (x - ai.toDouble())
    }
    if ((m00 / m10).toInt().also { number[0] = it } != 0) {
        number[1] = (m00 - number[0] * m10).toInt()
        number[0] *= sign
        number[2] = m10.toInt()
        return number
    }
    number[0] = 0
    number[1] = (sign * m00).toInt()
    number[2] = m10.toInt()
    return number
}

fun formatRational(value: Double, fullRational: Boolean = true): String {
    val number: IntArray = decimalToRational(value)
    if (number[1] == 0) {
        return number[0].toString()
    }
    if (number[0] == 0) {
        return number[1].toString() + "/" + number[2].toString()
    }
    return if (fullRational) {
        ((number[1] + abs(number[0]) * number[2]) * if (number[0] >= 0) 1 else -1).toString() + "/" + number[2].toString()
    } else number[0].toString() + " " + number[1].toString() + "/" + number[2].toString()
}

fun formatRational(value: Vector, fullRational: Boolean = true): String {
    val str = StringBuilder("{ ")
    for (i in 0 until value.size - 1) {
        str.append(formatRational(value[i], fullRational))
        str.append(", ")
    }
    str.append(formatRational(value[value.size - 1], fullRational))
    str.append(" }")
    return str.toString()
}

fun Vector.toRationalStr() = formatRational(this)