package com.ssau.lib

class Logger {
    companion object {
        @JvmStatic
        fun log(string: String) {
            println(string)
        }
        @JvmStatic
        fun logIters(tag: String, iters: Int) {
            log("$tag iterations number: $iters")
        }
    }
}