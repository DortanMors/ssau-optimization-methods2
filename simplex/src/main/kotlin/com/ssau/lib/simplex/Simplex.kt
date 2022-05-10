package com.ssau.lib.simplex

import org.jetbrains.kotlinx.multik.ndarray.data.get
import org.jetbrains.kotlinx.multik.ndarray.data.set
import org.jetbrains.kotlinx.multik.ndarray.operations.*
import kotlin.math.abs

class Simplex(initialBoundsMatrix: List<List<Double>>, initialPricesVector: List<Double>, initialInequalities: List<Sign>, initialBoundsVector: List<Double>) {

    private var boundsVector: Vector = vector(initialBoundsVector)
    private var pricesVector: Vector = vector(initialPricesVector)
    private var boundsMatrix: Matrix = matrix(initialBoundsMatrix)
    private var extremumMode: Extremum = Extremum.Min
    private lateinit var naturalArgsIds: MutableList<Int>
    private lateinit var basisArgs: MutableList<Int>
    private lateinit var fModArgs: MutableList<Int>
    private lateinit var artificialArgsIds: MutableList<Int>

    private lateinit var simplexTable: Matrix
    private var inequalities = initialInequalities.toMutableList()

    fun findMax() = apply { extremumMode = Extremum.Max }

    fun findMin() = apply { extremumMode = Extremum.Min }

    fun withBoundsMatrix(matrix: List<List<Double>>) = apply { boundsMatrix = matrix(matrix) }

    fun withBoundsVector(vector: List<Double>) = apply { boundsVector = vector(vector) }

    fun withPrices(vector: List<Double>) = apply { pricesVector = vector(vector) }

    fun withInequalities(list: List<Sign>) = apply { inequalities = list.toMutableList() }


    init {
        if (boundsVector.size != inequalities.size) {
            throw IllegalArgumentException("Error simplex creation :: b.size() != inequalities.size()")
        }
        if (boundsMatrix.rows() != inequalities.size) {
            throw IllegalArgumentException("Error simplex creation :: A.rows_number() != inequalities.size()")
        }
        if (boundsMatrix.cols() != pricesVector.size) {
            throw IllegalArgumentException("Error simplex creation :: A.cols_number() != price coefficients.size()")
        }
    }

    fun solve(): Vector? {
        println("Simplex finding $extremumMode")

        buildSimplexTable()

        var aik: Double

        var mainRow: Int

        var mainCol: Int

        println("Start simplex table:\n$this")

        if (excludeModArgs()) {
            // второй этап, если задача должна решаться двух проходным(двух этапным) алгоритмом
            println("Simplex table after args exclude:\n$this")
        }

        while (!isPlanOptimal()) {
            mainCol = getMainCol()
            if (mainCol == -1) {
                break
            }
            mainRow = getMainRow(mainCol)
            if (mainRow == -1) {
                // Невозможность определить ведущую строку свидейтельствует о том, что обрасть поиска неограничена
                println("Unable to get main row. Simplex is probably boundless...")
                return null
            }
            basisArgs[mainRow] = mainCol
            aik = simplexTable[mainRow][mainCol]
            simplexTable[mainRow] = simplexTable[mainRow] * (1.0 / aik)
            for (i in 0 until simplexTable.rows()) { // TODO rows size
                if (i == mainRow) {
                    continue
                }
                simplexTable[i] -= simplexTable[i][mainCol] * simplexTable[mainRow] // TODO ??
            }
            println("a_main{" + (mainRow + 1).toString() + ", " + (mainCol + 1).toString() + "} = " + formatRational(aik))
            println(this)
            println("current solution: " + formatRational(currentSimplexSolution()))
        }
        if (validateSolution()) {
            return currentSimplexSolution(true).also { solution ->
                println("solution: " + solution.toRationalStr())
            }
            // формирование ответа
        }
        println("Simplex is unresolvable")
        // значение целевой функции не равно ее значению от найденного плана
        return null
    }

    private fun buildSimplexTable() {
        simplexTable = boundsMatrix.copy()
        naturalArgsIds = ArrayList()
        basisArgs = ArrayList()
        artificialArgsIds = ArrayList()
        fModArgs = ArrayList()

        // Если среди вектора b есть отрицательные значения, то соответствующие строки
        // матрицы ограничений умножаем на мину один и меняем знак сравнения
        for (row in 0 until simplexTable.rows()) { // rows
            if (boundsVector[row] >= 0) {
                continue
            }
            inequalities[row] = if (inequalities[row] == Sign.Less) Sign.Greater else Sign.Less
            boundsVector[row] *= -1.0
            simplexTable[row] = simplexTable[row] * -1.0 // TODO getRow
        }
        for (i in 0 until pricesVector.size) {
            naturalArgsIds.add(i)
        }

        // построение искуственного базиса
        var basisArgsInfo: IntArray
        for (inequalityId in inequalities.indices) {
            basisArgsInfo = buildVirtualBasisCol(inequalityId, inequalities[inequalityId])
            naturalArgsIds.add(basisArgsInfo[0])
            if (basisArgsInfo[1] != -1) {
                basisArgs.add(basisArgsInfo[1])
                fModArgs.add(basisArgsInfo[1])
                artificialArgsIds.add(basisArgsInfo[1])
                continue
            }
            basisArgs.add(basisArgsInfo[0])
        }

        // добавим столбец ограницений
        simplexTable = simplexTable.withColumn(boundsVector)

        // Построение симплекс разностей
        val sDeltas: Vector = vector(simplexTable.cols())
        if (extremumMode == Extremum.Max) {
            for (j in sDeltas.indices) {
                sDeltas[j] = if (j < pricesVector.size) -pricesVector[j] else 0.0
            }
        } else {
            for (j in sDeltas.indices) {
                sDeltas[j] = if (j < pricesVector.size) pricesVector[j] else 0.0
            }
        }
        simplexTable = simplexTable.withRow(sDeltas) // TODO add row

        // Если целевая функуция не была модифицирована
        if (!isTargetFuncModified()) {
            return
        }

        // Если всё же была...
        val sDeltasAdd: Vector = vector(simplexTable.cols())
        for (j in fModArgs.indices) {
            sDeltasAdd[fModArgs[j]] = 1.0
        }
        simplexTable = simplexTable.withRow(sDeltasAdd)
    }

    private fun naturalArgsN() = pricesVector.size

    private enum class Extremum { Min, Max }

    override fun toString(): String {
        if (simplexTable.rows() == 0) {
            return ""
        }
        val sb = StringBuilder()
        var i = 0
        sb.append(String.format("%-6s", " "))
        while (i < simplexTable.cols() - 1) {
            sb.append(String.format("|%-12s", " x " + (i + 1).toString()))
            i++
        }
        sb.append(String.format("|%-12s", " b"))
        sb.append("\n")
        var nRow = -1
        for (row in 0 until simplexTable.shape[0]) { // TODO: rows
            nRow++
            if (isTargetFuncModified()) {
                when (nRow) {
                    simplexTable.rows() - 2 -> {
                        sb.append(String.format("%-6s", " d0"))
                    }
                    simplexTable.rows() - 1 -> {
                        sb.append(String.format("%-6s", " d1"))
                    }
                    else -> {
                        sb.append(String.format("%-6s", " x " + (basisArgs[nRow] + 1).toString()))
                    }
                }
            } else {
                if (nRow == simplexTable.rows() - 1) {
                    sb.append(String.format("%-6s", " d"))
                } else {
                    sb.append(String.format("%-6s", " x " + (basisArgs[nRow] + 1).toString()))
                }
            }
            for (col in simplexTable[row]) {
                if (col >= 0) {
                    sb.append(java.lang.String.format("|%-12s", " " + formatRational(col)))
                    continue
                }
                sb.append(java.lang.String.format("|%-12s", formatRational(col)))
            }
            sb.append("\n")
        }
        sb.append("\n")
        return sb.toString()
    }

    /**
     * Проверяет оптимальность текущего опорного плана.
     * Исследуется положительность симплекс-разностей в последней строке СТ в диапазоне от 1:n-1.
     * Если целевая функция была модифицирована, то исследует две последних строки.
     * Если среди элементов от 1:n-1 в последней строке нет отрицательных, то проверяет
     * на неотрицательность только те элементы предпоследней строки, которые не являются искусственными.
     * @return
     */
    private fun isPlanOptimal(): Boolean {
        // Проверяем значения последней строки сиплекс-разностей
        // на положительность. Если все положительны, то план оптимален
        val row = simplexTable[simplexTable.rows() - 1]
        var opt = row.toList().none { it < 0 }

        // если мы модифицировали целевую функцию, то среди списка естественнхых
        // агументов проверям на положительнность предпоследнюю строку симплекс-разностей
        if (isTargetFuncModified()) {
            if (!opt) {
                return false
            }
            val rowModified = simplexTable[simplexTable.rows() - 2]


            for (id in naturalArgsIds) {
                if (rowModified[id] < 0) {
                    opt = false
                    break
                }
            }
        }
        return opt
    }

    private fun isTargetFuncModified() = fModArgs.size != 0

    /**
     * строит виртуальный базисный вектор
     * @param inequalityId
     * @param inequalitySign
     * @return
     */
    private fun buildVirtualBasisCol(inequalityId: Int, inequalitySign: Sign): IntArray =
        when (inequalitySign) {
            Sign.Less,
            Sign.Equal -> {
                val column = (0 until simplexTable.rows()).map { row ->
                    if (row == inequalityId) 1.0 else 0.0
                }
                simplexTable = simplexTable.withColumn(column)
                val cols = if (inequalitySign == Sign.Equal) simplexTable.cols() - 1 else -1
                intArrayOf(simplexTable.cols() - 1, cols)
            }
            Sign.Greater -> {
                val col1 = ArrayList<Double>()
                val col2 = ArrayList<Double>()
                for (row in 0 until simplexTable.rows()) {
                    if (row == inequalityId) {
                        col1 += (-1.0)
                        col2 += (1.0)
                    } else {
                        col1 += (0.0)
                        col2 += (0.0)
                    }
                }
                simplexTable = simplexTable.withColumn(col1)
                simplexTable = simplexTable.withColumn(col2)
                intArrayOf(simplexTable.cols() - 2, simplexTable.cols() - 1)
            }
        }


    private fun validateSolution(): Boolean {
        var `val` = 0.0
        val nRows = if (isTargetFuncModified()) simplexTable.rows() - 2 else simplexTable.rows() - 1
        val nCols = simplexTable.cols() - 1
        for (i in basisArgs.indices) {
            if (basisArgs[i] < naturalArgsN()) {
                `val` += simplexTable[i, nCols] * pricesVector[basisArgs[i]]
            }
        }
        if (extremumMode == Extremum.Max) {
            if (abs(`val` - simplexTable[nRows, nCols]) < 1e-5) {
                return if (isTargetFuncModified()) {
                    abs(simplexTable[simplexTable.rows() - 1, simplexTable.cols() - 1]) < 1e-5
                } else true
            }
        }
        return if (abs(`val` + simplexTable[nRows, nCols]) < 1e-5) {
            if (isTargetFuncModified()) {
                abs(simplexTable[simplexTable.rows() - 1, simplexTable.cols() - 1]) < 1e-5
            } else true
        } else false
    }

    private fun excludeModArgs(): Boolean {
        if (!isTargetFuncModified()) {
            return false
        }
        val lastRowId = simplexTable.rows() - 1
        for (i in fModArgs.indices) {
            for (row in 0 until simplexTable.rows()) {
                if (simplexTable[row][fModArgs[i]] != 0.0) {
                    val arg = simplexTable[lastRowId, fModArgs[i]] / simplexTable[row, fModArgs[i]]
                    simplexTable[lastRowId] -= arg * simplexTable[row]
                    break
                }
            }
        }
        return true
    }

    /**
     * Определяет ведущий столбец. Среди элементов строки симплекс-разностей ищет максимальны по модулю
     * отрицательный элемент. Если целевая функция была модифицирована и среди последней строки нет отрицательных
     * элементов, то посик таковых будет продолжен среди только тех элементов предпоследней строки, которые не
     * являются искусственными.
     * @return
     */
    private fun getMainCol(): Int {
        val row = simplexTable[simplexTable.rows() - 1]
        var delta = 0.0
        var index = -1
        for (i in 0 until row.size - 1) {
            if (row[i] >= delta) {
                continue
            }
            delta = row[i]
            index = i
        }
        if (isTargetFuncModified() && index == -1) {
            val rowAdd = simplexTable[simplexTable.rows() - 2]
            for (id in naturalArgsIds) {
                if (rowAdd[id] >= delta) {
                    continue
                }
                delta = rowAdd[id]
                index = id
            }
        }
        return index
    }

    /**
     * Определяет ведущую строку
     * @param simplexCol ведущий столбец
     * @return
     */
    private fun getMainRow(simplexCol: Int): Int {
        var delta = 1e12
        var index = -1
        var aik: Double
        val bIndex = simplexTable.cols() - 1
        val rowsN = if (isTargetFuncModified()) simplexTable.rows() - 2 else simplexTable.rows() - 1
        for (i in 0 until rowsN) {
            aik = simplexTable[i, simplexCol]
            if (aik < 0) {
                continue
            }
            if (simplexTable[i, bIndex] / aik > delta) {
                continue
            }
            delta = simplexTable[i, bIndex] / aik
            index = i
        }
        return index
    }

    private fun currentSimplexSolution(onlyNaturalArgs: Boolean = false): Vector {
        val solution: Vector = vector(if (onlyNaturalArgs) naturalArgsN() else simplexTable.cols() - 1)
        for (i in 0 until basisArgs.size) {
            if (basisArgs[i] >= solution.size) {
                continue
            }
            solution[basisArgs[i]] = simplexTable[i, simplexTable.cols() - 1]
        }
        return solution
    }

}
