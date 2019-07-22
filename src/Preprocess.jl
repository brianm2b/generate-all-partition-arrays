# Preprocess.jl
"""
    isInputMatrixProper()

Determine whether or not the input matrix can fulfill the minimum non-musical requirements of an all-partition array.
"""
function isInputMatrixProper(inputMatrix)
    print("\n", "Begin preprocessing...")
    numRows = size(inputMatrix, 1)
    numCols = size(inputMatrix, 2)
    if isless(numCols, numRows)
        throw(DimensionMismatch("input matrix must have fewer rows than columns")) 
    end
    if numRows > 12 || numRows < 4
        throw(DimensionMismatch("input matrix must have more than 3 rows and fewer than 13"))
    end
    if numCols > 96 || numCols < 72
        throw(DimensionMismatch("input matrix must have more than 71 columns and fewer than 97"))
    end
    aggregate = sort(unique(inputMatrix))
    numPcs = length(aggregate)
    if numPcs != 12
        throw(ErrorException("number of distinct elements in input matrix must be equal to 12"))
    end
    if !isequal(aggregate, [0,1,2,3,4,5,6,7,8,9,10,11])
        throw(ErrorException("elements of input matrix must be in the range [0, 11]"))
    end
    counts = zeros(numPcs)
    for i in 1:numPcs
        counts[i] = count(isequal(aggregate[i]), inputMatrix)
    end
    if !all(counts[1] .== counts)
        throw(ErrorException("input matrix must contain an equal number of each element in the range [0, 11]"))
    end
    aggregate .+= 1 # add one to match matrix adjustment done later
    numElements::Int16 = length(inputMatrix)
    numPartitions = findNumOfPartitions(numRows, numPcs)
    numOverlaps::Int16 = (numPartitions * numPcs) - (numRows * numCols)
    if isless((convert(Int16, numPartitions) * convert(Int16, numPcs)), numElements) # FIX
        throw(ErrorException("input matrix must contain a number of elements equal to or fewer than the number required of its partitions"))
    end
    return APArray(numRows, numCols, numElements, numOverlaps, numPartitions, aggregate) 
end

"""
    findNumOfPartitions()

Find the number of required partitions given the number of rows in the input matrix and the size of the aggregate.
"""
function findNumOfPartitions(numRows, numPcs)
    numPartitions = zero(eltype(numRows))
    for i in 1:numRows
        numPartitions += length(collect(partitions(numPcs, i)))
    end
    return numPartitions
end

"""
    requestUserInputToSelectHeuristics()

Halts the continuation of the program until the user provides appropriate input to either continue or stop the program.
If the user opts to continue, they can select whether or not to use searhc heuristics. If they choose to use heuristics,
they can select one of two options.
"""
function requestUserInputToSelectHeuristics(arraySpecs, reservedPartitions, allCompositions)
    userChoiceOffset = zero(Int8)
    userChoiceHeuristic = Array{Int8,1}()
    userChoiceReserves = zero(Int8)
    println("The generateAllPartitionArray() algorithm will search through the input matrix and attempt to generate an all-partition array with the following specifications:","\n")
    printInputMatrixSummary(arraySpecs, reservedPartitions, allCompositions) # print details of input matrix and its an all-partition array
    # ask for user input to begin preprocessing and searching for a solution
    userInput = ""
    while (userInput == "")
        println("")
        userInput = getUserInput("This process can include the use of heuristics that work globally to help guide the search. To select from these heuristics and begin searching, press 'y' then RETURN. To terminate the program now, press 'n' then RETURN.")    
        #userInput = parse(AbstractString, getUserInput("This process can include several heuristics that work globally to help guide the search. To select from these heuristcs and begin searching, press 'y' then RETURN. To terminate the program now, press 'n' then RETURN."))    
        if isequal(userInput, "n")
            return exit()
        elseif !isequal(userInput, "y")
            println("Incorrect key entered", "\n") # use ErrorException?
            userInput = ""
        end
    end
    # ask for user input for selecting a heuristic or not
    userInput = ""
    while (userInput == "")
        println("")
        # ask for user input for selecting offset
        userInput = getUserInput("Input ONE of the following numbers: '0' for no search heuristic (not recommended), or '1' for use of a search heuristic. Afterwards, press RETURN.")   
        userChoiceHeuristic = [parse(Int8,i) for i in userInput]
        if length(userChoiceHeuristic) == 0
            println("Incorrect key entered", "\n")
            userInput = ""
        end
    end
    if userChoiceHeuristic[1] != 0
        println("\n", "Select the options below if you are using any of the following input matrix types:","\n")
        println("Babbitt(6,96):","\t","1 search heuristic and ", 13, " (i.e., maximum)", " reserved partitions","\n")
        println("Smalley(6,96):","\t","0 search heuristic and ", 22, " (i.e., maximum)", " reserved partitions","\n")
        println("Babbitt(4,96):","\t","1 search heuristic and ", 1, " reserved partitions","\n")

        # ask for user input for selecting the search heuristic
        userInput = ""
        while (userInput == "")
            println("")
            userInput = getUserInput("Input ONE of the following numbers for judging candidates with respect to an 'ideal position': '0' for only the magnitude of difference, or '1' for the individual magnitude and direction of difference. Afterwards, press RETURN.")   
            userChoiceOffset = parse(Int8, userInput)
            if length(userChoiceOffset) != 1
                println("Input must be a single number", "\n")
                userInput = ""
            end
            if length(userChoiceOffset) == 0
                println("Incorrect key entered", "\n")
                userInput = ""
            end
        end
    end
    # ask for user input for selecting reserved partitions 
    userInput = ""
    while (userInput == "")
        println("")
        userInput = getUserInput("Input the number of final region partitions––from 0 to the number of possible reserved partitions (above) to reserve in the search. Afterwards, press RETURN.")   
        userChoiceReserves = parse(Int8, userInput)
        if userChoiceReserves > length(reservedPartitions) || userChoiceReserves < 0
            println("Input must be a number between 0 and ", length(reservedPartitions),"\n")
            userInput = ""
        end
    end
    return Heuristic(userChoiceOffset, userChoiceHeuristic, userChoiceReserves)
end

"""
    printInputMatrixSummary()

"""
function printInputMatrixSummary(arraySpecs, reservedPartitions, allCompositions) # print more info
    println("Input matrix size:", "\t", "\t", "\t", "\t", arraySpecs.numRows, " x ", arraySpecs.numCols)
    println("Aggregate size:", "\t", "\t", "\t", "\t", "\t", length(arraySpecs.aggregate))
    println("Num. of required partitions:", "\t", "\t", "\t", arraySpecs.numPartitions, "\t", "(from ", size(allCompositions, 1), " compositions)")
    println("Num. of possible reserved partitions:", "\t", "\t", length(reservedPartitions))
    println("Num. of require overlaps:", "\t", "\t", "\t", arraySpecs.numOverlaps)
    #println("Summary: Input matrix is proper.", "\n")
    return nothing
end

"""
    getUserInput()

Prompt user for input and return user input as a string.
"""
function getUserInput(prompt::AbstractString="")
    println(prompt * " ")
    convert(String, chomp(readline()))
end

"""
    computeListOfAllCompositions()

Computes a list of all compositions needed for the given input matrix.
Also finds the sum of the largest and smallest parts of all partitions and creates a helper dictionary variable for searching all compositions.
"""
function computeListOfAllCompositions(arraySpecs)
    lexCompositions = computeLexicalOrderCompostitions(length(arraySpecs.aggregate), arraySpecs.numRows)
    allCompositions, dictPartComps, minMaxParts, probabilities = groupCompositionsByPartition(lexCompositions, arraySpecs.numRows, arraySpecs.numPartitions)
    return allCompositions, dictPartComps, minMaxParts, probabilities
end

"""
    computeLexicalOrderCompostitions()

Computes a lexicographically ordered list of all compositions with n parts or fewer that sum to k, where n is the number 
of input matrix rows and k is the number of pitch classes in an aggregate.
"""
function computeLexicalOrderCompostitions(numPcs, numRows)
    # k: sum, n: number of non-negative integers (i.e., n sum k)
    m = binomial(numPcs + numRows - 1, numRows - 1)
    d1 = zeros(eltype(m), m, 1)
    d2 = collect(combinations(collect(1:numPcs + numRows - 1), numRows - 1))
    d2 = transpose(hcat(d2...)) #
    d3 = fill!(ones(eltype(m), m, 1), numPcs + numRows)
    dividers = hcat(d1, d2, d3)
    return convert(Array{Int8,2}, diff(dividers, dims=2) .- 1)
end

"""
    groupCompositionsByPartition()

Groups compositions by partition ordered lexicographically. Returns this list in addition to some helper variables.
"""
function groupCompositionsByPartition(lexCompositions, numRows, numPartitions)
    allCompositions = zeros(eltype(lexCompositions), size(lexCompositions, 1), size(lexCompositions, 2)+1)
    dictPartComps = Dict{eltype(lexCompositions), UnitRange}()
    minMaxParts = zeros(Int16, numPartitions, 2)
    uniqueMatrix = unique(sort(lexCompositions, dims=2), dims=1)
    probabilities = zeros(Float16, numPartitions)
    start = 1
    for i in 1:numPartitions
        testComp = uniqueMatrix[i, :]
        indices = Array{Int,1}()
        for j in 1:size(lexCompositions, 1)
            if testComp == sort(lexCompositions[j, :])
                push!(indices, j)
            end
        end
        allCompositions[start:start+length(indices)-1, 1:numRows] = lexCompositions[indices, 1:numRows] # store composition
        allCompositions[start:start+length(indices)-1, end] .= i # store partition number
        dictPartComps[i] = start:start+length(indices)-1 # create dictionary entry of key: partition number, value: composition index range
        minMaxParts[i, 1] = first(testComp) # store smallest part
        minMaxParts[i, 2] = last(testComp) # store largest part
        # Compute probabilities here
        #probabilities[i] = length(indices) / size(allCompositions, 1) # number of compositions
        probabilities[i] = length(testComp) - count(iszero, testComp)  # number of zeros in composition
        start = start + length(indices)
    end
    # normalize probabilities
    probabilities = (probabilities .- minimum(probabilities))/(maximum(probabilities) - minimum(probabilities))
    return allCompositions, dictPartComps, minMaxParts, probabilities
end


"""
    computeListOfAllCombinationsOfOverlaps()

Find all possible overlaps for a composition with at most a number of parts equal to the number of rows 
in the input matrix. The number of possible overlaps is equal to all combinations of the number of matrix rows.
"""
function computeListOfAllCombinationsOfOverlaps(numRows)
    indices = collect(Array{Int16,1}, combinations(1:numRows))
    allOverlaps = falses(length(indices) + 1, numRows)
    # rowsOverlapsByNum = [1, 7, 22, 42, 57, 63] # e.g., greater than x[1] is 1 or more overlaps 
    rowsOverlapsByNum = Array{Int16,1}()
    setSize = 0
    for i::Int16 in 2:length(indices) + 1
        allOverlaps[i, indices[i - 1]] .= true
        if length(indices[i - 1]) != setSize
            push!(rowsOverlapsByNum, i - one(eltype(rowsOverlapsByNum))) 
            setSize = length(indices[i - 1])
        end
    end
    overlapIndicesByRow = findOverlapRowIndices(allOverlaps)
    return allOverlaps, rowsOverlapsByNum, overlapIndicesByRow
end

"""
    findOverlapRowIndices()

Find all rows of overlaps that share one overlap and return a list of the indices of these rows in 
the set of all combinations of overlaps (one list for each input matrix row).
"""
function findOverlapRowIndices(allOverlaps)
    # group row indices of overlaps by shared rows (12-row input matrix: max 4095)
    overlapIndicesByRow = zeros(Int16, Int(size(allOverlaps, 1) / 2), size(allOverlaps, 2))
    for i in 1:size(allOverlaps, 2)
        overlapIndicesByRow[:, i] = findall(allOverlaps[:, i])
    end
    return overlapIndicesByRow
end

"""
    findPartitionsToReserveForLast()

Finds all candidate partitions at the righthand side of the imput matrix (from right to left) corresponding 
to the possible partitions for the final position in the primary left-to-right search. These partitions can
be excluded from the primary search so that at least one remains available for selection for the final position.
"""
function findPartitionsToReserveForLast(inputMatrix, allCompositions, allOverlaps, overlapIndicesByRow, arraySpecs, minMaxParts, rowsOverlapsByNum, probabilities)
    matrix = reverse(inputMatrix, dims=2) # work backwards from matrix
    matrix .+= one(eltype(inputMatrix))
    position = ones(eltype(inputMatrix), arraySpecs.numRows)   # matrix columns for each row
    cnt = one(eltype(inputMatrix))
    usedPartitions = zeros(Int, arraySpecs.numPartitions)
    sumUnUsedParts = computeSumOfUnUsedParts(usedPartitions, minMaxParts, cnt)
    idealRegion::Float16 = round((arraySpecs.numRows * (arraySpecs.numCols + 1)) / arraySpecs.numPartitions, digits=3)
    chosenHeuristics = Heuristic(0, Array{Int8,1}([3]), 0)

    # find final position compositions
    # two dummy variables (x, y)
    candidateList, x, y = findCandidateCompositionsAtPosition(matrix, position, cnt, allCompositions, allOverlaps, overlapIndicesByRow, arraySpecs, minMaxParts, sumUnUsedParts, rowsOverlapsByNum, chosenHeuristics, idealRegion, probabilities)
    
    # find partitions from the final position compositions
    reservedLocs = falses(arraySpecs.numPartitions)
    for i in 1:length(candidateList)
        reservedLocs[candidateList[i].partition] = true
    end
    reservedPartitions::Array{eltype(inputMatrix),1} = 1:arraySpecs.numPartitions
    reservedPartitions = reservedPartitions[reservedLocs]
    
    #= Request user input to select heuristics used in searching =#
    print("[complete]")
    println("\n")
    chosenHeuristics = requestUserInputToSelectHeuristics(arraySpecs, reservedPartitions, allCompositions) # user-selected heuristics

    return reservedPartitions, chosenHeuristics
end

#= End of file =#