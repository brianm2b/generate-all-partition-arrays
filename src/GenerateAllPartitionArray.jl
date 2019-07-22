# GenerateAllPartitionArray.jl
module GenerateAllPartitionArray

    using Combinatorics, DelimitedFiles, Statistics

    export 
        # preprocessing functions
        isInputMatrixProper, 
        findNumOfPartitions, 
        requestUserInputToSelectHeuristics, 
        printInputMatrixSummary, 
        getUserInput, 
        computeListOfAllCompositions, 
        computeLexicalOrderCompostitions, 
        groupCompositionsByPartition, 
        computeListOfAllCombinationsOfOverlaps, 
        findOverlapRowIndices, 
        findPartitionsToReserveForLast,

        # core functions
        generateAllPartitionArray,
        searchMatrixForCoveringOfCandidates,
        makeInputMatrixSuitableForSearch,
        findUnUsedCompositions, 
        findUnUsedPartitions,
        findSumOfUnUsedParts, 
        findCandidateCompositionsAtPosition, 
        isCompositionPossible, 
        isRegionPossible, 
        findMissingRegionPcs, 
        findDuplicateRegionPcs, 
        findPossibleOverlapsForRegion, 
        findCandidateOverlapsForRegion, 
        selectHeuristicsAndCompute,
        judgeQualityOfCandidate,
        sortCandidatesByHeuristicScore,
        success, 
        backtrack, 
        printSolution

        # files containing preprocessing and core functions
        include("./Preprocess.jl")
        include("./Core.jl")


#= Read in input matrix from file =#
const inputMatrix = readdlm("Babbitt_6_96.txt", '\t', Int8, '\n')   # as found e.g., in Babbitt's Arie da Capo
#const inputMatrix = readdlm("Smalley_6_96.txt", '\t', Int8, '\n')  # as found e.g., in Babbitt's Sheer Pluck (constructed by David Smalley)
#const inputMatrix = readdlm("Babbitt_4_96.txt", '\t', Int8, '\n')  # as found e.g., in Babbitt's My Ends are My Beginings


#= Define primary types =#
"""
    Immutable composite type called Candidate

A left-to-right sequence of K Candidates are needed to form a cover (i.e., a solution and complete all-partition array).        
"""
struct Candidate
    partition::Int8
    composition::Array{Int8,1}
    overlaps::BitArray{2}
end

"""
    Mutable composite type called Region
"""
mutable struct Region
    pcs::Array{Int8,1}
    inPushPcs::Array{Int8,1} # indexed by overlapsAllowedHere
    outPushPcs::Array{Int8,1} # indexed by overlapsAllowedHere
    missingPcs::Array{Int8,1} #
    duplicatePcs::Array{Int8,1} #
    # inner constructor of default values
    # Region() = new(zeros(Int8,12), zeros(Int8,6), zeros(Int8,6), Array{Int8,1}(), Array{Int8,1}()) # remove?
end

"""
    Immutable composite type called APArray    
"""
struct APArray
    numRows::Int8
    numCols::Int8
    numElements::Int16
    numOverlaps::Int16
    numPartitions::Int8
    aggregate::Array{Int8,1}
end

"""
    Immutable composite type called Heuristic     
"""
struct Heuristic
    offset::Int8
    heuristics::Array{Int8,1}
    reserves::Int8
end

"""
    Mutable composite type called Selected
"""
mutable struct Selected
    candidates::Array{Int32,1}
    overlaps::Array{Int16,1}
end

#= Run main algorithm =#
selected, allCandidateLists, cnt, position = generateAllPartitionArray(inputMatrix)

end
#= End of module =#
