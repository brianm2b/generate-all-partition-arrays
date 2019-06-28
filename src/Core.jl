# Core.jl
"""
    generateAllPartitionArray()

Generates a musical structure created by Milton Babbitt (1916–2011) known as an all-partition array. An all-parititon array 
is an N-by-M regular matrix, A, containing (N * M) / 12 integers each of [0, 11] that can be "covered" by the K distinct 
partitions of 12 into N parts or fewer. Each of these K partitions contains 12 distinct matrix elements [0, 11] whose N 
or fewer parts each consist of consecutive row locations in A. Parts from each of these K partitions are "packed" or ordered 
sequentially in A from left to right. The number of elements in these partitions, K * 12, required to cover A exceeds the 
number of elements in A. For this reason, an all-partition array will have (K * 12) - (N * M) overlapping locations in A 
between contiguous partitions. Such partitions may overlap in at most one row location in A for each of their contiguous 
parts. This row location in A must be the right-most element of some part in a partition k and the left-most element of 
a part in the same row in A of some contiguous partition greater than k (i.e., to its right).
"""
function generateAllPartitionArray(inputMatrix)
    #= Begin preprocessing =#
    arraySpecs = isInputMatrixProper(inputMatrix) # check if input matrix fulfills necessary requirements of an all-partition array
    allCompositions, dictPartComps, minMaxParts, probabilities = computeListOfAllCompositions(arraySpecs)
    allOverlaps, rowsOverlapsByNum, overlapIndicesByRow = computeListOfAllCombinationsOfOverlaps(arraySpecs.numRows)
    reservedPartitions, chosenHeuristics = findPartitionsToReserveForLast(inputMatrix, allCompositions, allOverlaps, overlapIndicesByRow, arraySpecs, minMaxParts, rowsOverlapsByNum, probabilities)

    #= Begin searching =#
    selected, allCandidateLists, cnt, position, numBacktracks = searchMatrixForCoveringOfCandidates(inputMatrix, allCompositions, reservedPartitions, allOverlaps, dictPartComps, rowsOverlapsByNum, overlapIndicesByRow, minMaxParts, arraySpecs, chosenHeuristics, probabilities) 
    printSearchResults(selected, allCandidateLists, cnt, inputMatrix, numBacktracks)

    return selected, allCandidateLists, cnt, position
end

#= Search matrix for candidate compostitions =#
"""
    searchMatrixForCoveringOfCandidates() A.K.A. the "BacktrackingBabbitt" algorithm


"""
function searchMatrixForCoveringOfCandidates(inputMatrix, allCompositions, reservedPartitions,
                                            allOverlaps, dictPartComps, rowsOverlapsByNum, 
                                            overlapIndicesByRow, minMaxParts, arraySpecs, chosenHeuristics,
                                            probabilities) 
    #= Primary search variables =#
    searchMatrix = makeInputMatrixSuitableForSearch(inputMatrix) # pad input matrix and add 1 to each element
    cnt = one(eltype(searchMatrix)) # count of the number of partitions used
    position = ones(eltype(searchMatrix), arraySpecs.numRows)  # row positions in matrix for a given cnt
    stopPosition = fill(arraySpecs.numCols + 1, arraySpecs.numRows) # row positions when solution found (check type)
    selected = Selected(zeros(Int32, arraySpecs.numPartitions), zeros(Int16, arraySpecs.numPartitions))
    usedPartitions = zeros(eltype(searchMatrix), arraySpecs.numPartitions) # partition used at each cnt 
    allCandidateLists = Array{Candidate}[[] for i in 1:arraySpecs.numPartitions] # list of candidates found at each cnt
    farthestReached = zero(eltype(searchMatrix))
    numBacktracks = zero(Int64)
    #= Heuristic variables =# 
    idealRegion::Float16 = round((arraySpecs.numRows * (arraySpecs.numCols)) / arraySpecs.numPartitions, digits=3)
    #idealRegion::Float16 = round((arraySpecs.numRows * (arraySpecs.numCols + 1)) / arraySpecs.numPartitions, digits=3)
    
    #= Begin backtracking search =# 
    println("\n", "Begin searching...", "\n")
    @time while cnt <= arraySpecs.numPartitions && isless(0, cnt)
        if isempty(allCandidateLists[cnt])
            unUsedCompositions, sumUnUsedParts = findUnUsedCompositions(allCompositions, dictPartComps, reservedPartitions, usedPartitions, cnt, minMaxParts, chosenHeuristics)
            allCandidateLists[cnt], selected.candidates[cnt], selected.overlaps[cnt] = findCandidateCompositionsAtPosition(searchMatrix, position, cnt, unUsedCompositions, allOverlaps, overlapIndicesByRow, arraySpecs, minMaxParts, sumUnUsedParts, rowsOverlapsByNum, chosenHeuristics, idealRegion, probabilities)
            if isempty(allCandidateLists[cnt]) # no list of candidates found at this cnt for current position
                # backtrack
                if isless(1, cnt)
                    #cnt, position, allCandidateLists[cnt + one(eltype(cnt))], usedPartitions[cnt + one(eltype(cnt))], selected.candidates[cnt + one(eltype(cnt))], selected.overlaps[cnt + one(eltype(cnt))] = backtrack(position, cnt, allCandidateLists[cnt - one(eltype(cnt))][selected.candidates[cnt - one(eltype(cnt))]], selected.overlaps[cnt - one(eltype(cnt))])
                    cnt, position, allCandidateLists, usedPartitions, selected = backtrack(cnt, position, allCandidateLists, usedPartitions, selected)
                else
                    cnt, allCandidateLists[1], usedPartitions[1], selected.candidates[1], selected.overlaps[1] = backtrack(cnt)
                end
                numBacktracks += 1
            else # list of candidates found at this cnt for current position
                cnt, position, usedPartitions[cnt - one(eltype(cnt))] = success(position, cnt, allCandidateLists[cnt][selected.candidates[cnt]], selected.overlaps[cnt])
                #println(selected.candidates)
            end
        else # backtracked to a previously computed list of candidates at this cnt
            selected.overlaps[cnt] += one(eltype(selected.overlaps)) # choose next overlaps for previously selected candidate (may or may not have candidate overlaps)
            if isless(size(allCandidateLists[cnt][selected.candidates[cnt]].overlaps, 1), selected.overlaps[cnt]) # previously selected candidate has no remaining overlaps to select
                selected.candidates[cnt] += one(eltype(selected.candidates)) # select next available candidate in this list of candidates
                if isless(size(allCandidateLists[cnt], 1), selected.candidates[cnt]) # list has no remaining candidates to select
                    # backtrack
                    if isless(1, cnt)
                        #cnt, position, allCandidateLists[cnt + one(eltype(cnt))], usedPartitions[cnt + one(eltype(cnt))], selected.candidates[cnt + one(eltype(cnt))], selected.overlaps[cnt + one(eltype(cnt))] = backtrack(position, cnt, allCandidateLists[cnt - one(eltype(cnt))][selected.candidates[cnt - one(eltype(cnt))]], selected.overlaps[cnt - one(eltype(cnt))])
                        cnt, position, allCandidateLists, usedPartitions, selected = backtrack(cnt, position, allCandidateLists, usedPartitions, selected)
                    else
                        cnt, allCandidateLists[1], usedPartitions[1], selected.candidates[1], selected.overlaps[1] = backtrack(cnt)
                    end
                    numBacktracks += 1
                else # list has at least one remaining candidate (selected above)
                    selected.overlaps[cnt] = one(eltype(selected.overlaps)) # select first overlaps of newly selected candidate
                    cnt, position, usedPartitions[cnt - one(eltype(cnt))] = success(position, cnt, allCandidateLists[cnt][selected.candidates[cnt]], selected.overlaps[cnt])
                    #println(selected.candidates)
                end
            else # previously selected candidate has at least one remaining overlaps to select
                cnt, position, usedPartitions[cnt - one(eltype(cnt))] = success(position, cnt, allCandidateLists[cnt][selected.candidates[cnt]], selected.overlaps[cnt])
                #println(selected.candidates)
            end
        end
        farthestReached = printProgress(selected, farthestReached, cnt) # print furthest position reached 
        if isequal(cnt, arraySpecs.numPartitions + 1) # all partitions used
            if !isequal(position, stopPosition) # not a solution: matrix contains uncovered elements
                # special backtrack (returns 4 dummy variables in addition to cnt and position)
                cnt, position, w, x, y, z,  = backtrack(position, cnt, allCandidateLists[cnt - one(eltype(cnt))][selected.candidates[cnt - one(eltype(cnt))]], selected.overlaps[cnt - one(eltype(cnt))])
                numBacktracks += 1
            end
        end
    end
    return selected, allCandidateLists, cnt, position, numBacktracks
end

#= Make input matrix suitable for search =#
"""
    makeInputMatrixSuitableForSearch()

"""
function makeInputMatrixSuitableForSearch(inputMatrix)
    searchMatrix = copy(inputMatrix)
    searchMatrix .+= one(eltype(searchMatrix)) # no zeros allows for elements as indices in 1-based indexing
    searchMatrix = hcat(searchMatrix, fill(13, size(searchMatrix, 1), 13)) # consider rotating for efficient slicing
    return searchMatrix
end


#= Determine unused compositions =#
"""
    findUnUsedCompositions()

Compute list of all compositions belonging to the unused partitions in the left-to-right 
partial sequence of candidates chosen during the search process. Uses a dictionary.
"""
function findUnUsedCompositions(allCompositions, dictPartComps, reservedPartitions, usedPartitions, cnt, minMaxParts, chosenHeuristics)
    unUsedPartitions = findUnUsedPartitions(reservedPartitions, usedPartitions, cnt, chosenHeuristics)
    indices = falses(size(allCompositions, 1)) 
    for i in unUsedPartitions
        indices[dictPartComps[i]] .= true
    end
    return getindex(allCompositions, indices, :), computeSumOfUnUsedParts(usedPartitions, minMaxParts, cnt)
end

"""
    findUnUsedPartitions()

Determine which partitions do not appear in the left-to-right partial sequence of chosen candidates. 
Return a sorted list of K or fewer integers in the range, [1, K].
"""
function findUnUsedPartitions(reservedPartitions, usedPartitions, cnt, chosenHeuristics)
    indices = trues(length(usedPartitions))
    if isless(1, cnt) # partitions previously found (may or may not include any reserved partitions)
        indices[usedPartitions[1:cnt - one(eltype(cnt))]] .= false
    end
    # reserve r number of partitions until the end based on chosen value of r
    if cnt < (length(usedPartitions) - chosenHeuristics.reserves) + one(eltype(cnt)) && length(setdiff(reservedPartitions, usedPartitions[1:cnt - one(eltype(cnt))])) < chosenHeuristics.reserves + one(eltype(cnt))
       indices[reservedPartitions] .= false
    end
    return getindex(1:length(usedPartitions), indices)  # sorted
end

"""
    computeSumOfUnUsedParts()

Compute two integers: (1) the sum of the largest part from each of the unused partitions and 
(2) the sum of the smallest part from each of the unused partitions. In (2), the smallest part from 
each partition that is greater than zero has 1 subtracted from its value. This subtraction assumes 
the presence of an overlap.
"""
function computeSumOfUnUsedParts(usedPartitions, minMaxParts, cnt)
    indices = trues(length(usedPartitions))
    indices[usedPartitions[1:cnt - one(eltype(cnt))]] .= false  # include all unused and reserved partitions
    sumUnUsedSmallestParts = sum(minMaxParts[indices, 1]) - count(!iszero, minMaxParts[indices, 1]) # each non-zero part is assumed an overlap
    sumUnUsedLargestParts = sum(minMaxParts[indices, 2]) # no overlaps
    return vcat(sumUnUsedSmallestParts, sumUnUsedLargestParts)
end

#= Determine which of the unused compositions are candidates =#
"""
    findCandidateCompositionsAtPosition()

Determine which (if any) compositions belonging to the unused partitions at a given position k are 
candidates. A candidate is a composition which (1) does not make it impossible for the remaining parittions 
to fully cover the input matrix (or cover more), and (2) has a region that is either a complete aggregate or has missing and duplicate 
pcs that can be pushed to make it complete. 
"""
function findCandidateCompositionsAtPosition(searchMatrix, position, cnt, 
                                    unUsedCompositions, allOverlaps, overlapIndicesByRow, 
                                    arraySpecs, minMaxParts, sumUnUsedParts, rowsOverlapsByNum, 
                                    chosenHeuristics, idealRegion, probabilities)
    #= Primary variables =#
    currentRegion = Region(zeros(eltype(searchMatrix), length(arraySpecs.aggregate)), 
                            zeros(eltype(searchMatrix), length(position)), zeros(eltype(searchMatrix), 
                            length(position)), Array{eltype(searchMatrix),1}(), Array{eltype(searchMatrix),1}())
    candidateList = Array{Candidate,1}()
    candidateScores = Array{Array{Float16,1},1}()
    
    #= Heuristic variables =#
    # check type
    numUnusedPartElems = (arraySpecs.numPartitions * length(arraySpecs.aggregate)) - (length(arraySpecs.aggregate) * (cnt - 1))
    numUncoveredMatrixElems = (arraySpecs.numRows * arraySpecs.numCols) - sum(position) - length(position) # double check
    numUnusedOverlaps = arraySpecs.numOverlaps - (cnt * length(arraySpecs.aggregate)) - sum(position) - length(position)

    idealPos = findIdealPosition(cnt, position, chosenHeuristics, arraySpecs)

    #= Search for candidates =#
    for i in 1:size(unUsedCompositions, 1) #@code_warntype
        composition = getindex(unUsedCompositions, i, :) # includes partition number. inefficient row slice 
        possible, currentRegion, overlapsAllowedHere = isCompositionPossible(searchMatrix, position, composition, currentRegion, minMaxParts, sumUnUsedParts, arraySpecs.numCols)
        if possible # composition does not violate any constraints
            possible, currentRegion = isRegionPossible(currentRegion, composition, position, arraySpecs.aggregate, overlapsAllowedHere, numUnusedPartElems, arraySpecs.numCols) 
            if possible # region is either an aggregate or has missing and duplicate pcs that can be pushed to form an aggregate
                # all sets of possible overlaps that result in a complete region
                candidateOverlaps = findCandidateOverlapsForRegion(currentRegion, overlapsAllowedHere, allOverlaps, overlapIndicesByRow, rowsOverlapsByNum)
                if !isempty(candidateOverlaps) # candidate composition found
                    # heuristics for judging candidates and sorting their overlaps (if any)
                    score, candidateOverlaps = selectHeuristicsAndCompute(chosenHeuristics, position, idealPos, idealRegion, composition, candidateOverlaps, arraySpecs, probabilities[composition[end]])
                    push!(candidateScores, score)
                    push!(candidateList, Candidate(last(composition), composition[1:length(position)], candidateOverlaps)) # store candidate composition
                end
            end
            currentRegion.missingPcs = Array{eltype(searchMatrix),1}() # change this?
            currentRegion.duplicatePcs = Array{eltype(searchMatrix),1}()
        end
    end
    # check if no candidates found
    if isempty(candidateList)
        return candidateList, zero(Int32), zero(Int16)
    end
    # Break ties in rank using probs?
    return candidateList[sortCandidatesByHeuristicScore(candidateScores)], one(Int32), one(Int16) # sorted candidates by heuristic score
end

# 3 idealpositions: 1.66 for all, difference (+/-) from 1.66 for each, and avg. difference (+/-) from 1.66 for all
"""
    findIdealPosition()

Determines the "ideal" row positions of a sequence of candidates up to a given cnt. An ideal position can be one of three possibilities:
(1) (2) or (3)
"""
function findIdealPosition(cnt, position, chosenHeuristics, arraySpecs)
    # offsets is how far off the previous candidate was from ideal
    offsets = position - fill!(ones(Float16, length(position)), (cnt - 1 * round((arraySpecs.numCols) / arraySpecs.numPartitions, digits=3)))
    #if cnt == 1
    #    offsets = zeros(Float16,length(position))
    #end
    if isequal(chosenHeuristics.offset, 0) # constant ideal position (as used in BemmanMeredithJNMR2016, heuristic D in BemmanMeredithISMIR2019)
        # ideal position after having chosen a candidate at this cnt
        idealPos = fill!(ones(Float16, arraySpecs.numRows), (cnt * round((arraySpecs.numCols)/ arraySpecs.numPartitions, digits=3))) # 6-row: cnt * 1.66
    end
    if isequal(chosenHeuristics.offset, 1) # difference in magnitude and direction for each row (heuristic S in BemmanMeredithISMIR2019)
        idealPos = fill!(ones(Float16, arraySpecs.numRows), (cnt * round(((arraySpecs.numCols) / arraySpecs.numPartitions), digits=3))) - offsets #
    end
     #if isequal(chosenHeuristics.offset, 2) # average difference in magnitude and direction (i.e., same for all rows)
    #    idealPos = fill!(ones(Float16, arraySpecs.numRows), (cnt * round(((arraySpecs.numCols) / arraySpecs.numPartitions) - mean(offsets), digits=3))) #
    #end
    return idealPos
end


"""
    isCompositionPossible()

Determines whether or not a composition at a given position is possible. A composition is possible at a given position if 
(1) the remaining largest parts in unused partitions do not fail to cover the remaining pcs in any one matrix row and 
(2) the remaining smallest parts in unused partitions do not exceed the remaining pcs in any one matrix row. If a compostition
is possible, the function returns (1) the region of pcs it forms in the matrix, (2) the set of pcs lying to the left of each part,
(3) the set of pcs lying at the rightside edge of a part, and (4) the row location of where overlaps are allowed. 
"""
function isCompositionPossible(searchMatrix, position, composition, currentRegion, minMaxParts, sumUnUsedParts, numCols)
    overlapsAllowedHere = falses(length(position)) # assume no rows allowed overlaps
    satisfied, currentRegion, overlapsAllowedHere = searchMatrixRowsForRegion(searchMatrix, position, composition, currentRegion, overlapsAllowedHere, minMaxParts, sumUnUsedParts, numCols)
    satisfied || return false, currentRegion, overlapsAllowedHere # if not satisfied then return false
    return true, currentRegion, overlapsAllowedHere
end

"""
    searchMatrixRowsForRegion()

Finds several fields of an object of type Region() and the locations of allowed overlaps for a given composition, returning these if
the composition does not violate any constraints governing possible solutions.
"""
function searchMatrixRowsForRegion(searchMatrix, position, composition, currentRegion, overlapsAllowedHere, minMaxParts, sumUnUsedParts, numCols)
    prefix = one(eltype(composition))
    @inbounds for j::Int8 in 1:length(position) # consider switching to column-major order slicing
        if isless(0, composition[j]) # there is a part here
            if !compositionSatisfiesUnUsedPartitionsConstraint(position, composition, minMaxParts, sumUnUsedParts, numCols, j)
                return false, currentRegion, overlapsAllowedHere
            end
            currentRegion.pcs[prefix:prefix + composition[j] - one(eltype(composition))] = searchMatrix[j, position[j]:position[j] + composition[j] - one(eltype(composition))]
            if isless(1, position[j]) # overlap allowed in this matrix row
                setindex!(currentRegion.inPushPcs, searchMatrix[j, position[j] - one(eltype(position))], j)
                setindex!(currentRegion.outPushPcs, currentRegion.pcs[prefix + composition[j] - one(eltype(composition))], j)
                setindex!(overlapsAllowedHere, true, j)
            end
        end
        prefix += composition[j]
    end
    return true, currentRegion, overlapsAllowedHere
end

"""
    doesCompositionSatisfyRemainingPartitionsConstraint()

Determines whether or not the composition at a given position satisfies constraints governing possible solutions given the remaining number of 
partitions.
"""
function compositionSatisfiesUnUsedPartitionsConstraint(position, composition, minMaxParts, sumUnUsedParts, numCols, j)::Bool
    if isless(position[j] + composition[j] + last(sumUnUsedParts) - minMaxParts[composition[end],2], numCols + one(eltype(composition)))
        # not been enough parts (or too many overlaps) used in a single row
        return false
    end
    # could this fail if position has 1 (i.e., no overlaps allowed)?
    if isless(numCols + one(eltype(composition)), position[j] + composition[j] + first(sumUnUsedParts) - minMaxParts[composition[end],1] - one(eltype(composition))) 
        # not been enough overlaps (or too many parts) used in a single row
        return false
    end
    # tests the curent position NOT the composition
    #if isless(position[j] + last(sumUnUsedParts), numCols + one(eltype(composition)))
        # not been enough parts (or too many overlaps) used in a single row
    #    return false
    #end
    #if isless(last(sumUnUsedParts) - composition[j], numCols - position[j] + composition[j] - one(eltype(composition))) # subtract overlap?
        # not been enough parts (or too many overlaps) used in a single row
    #    return false
    #end
    return true
end

"""
    isRegionPossible()

Determine whether or not a region formed by a given composition at a given position is possible. A region is possible
if (1) there are at least as many possible overlaps as there are missing pcs from a region, (2) the number uncovered 
matrix elements does not exceed the number of unused partition elements, (3) it is not missing any pcs, (4) any missing
pcs are found in the set of pcs possible for "pushing" into, and (5) any duplicate pcs are found in the set of pcs
possible for "pushing" out.
"""
function isRegionPossible(currentRegion, composition, position, aggregate, overlapsAllowedHere, numUnusedPartElems, numCols)  
    missingPcHits = trues(length(aggregate) + 1) # hot-or-not occurences of 12 pcs plus possible 13th (matrix padded value "13")
    missingPcHits[currentRegion.pcs] .= false   # find locations of missing pcs
    # region is not possible if not enough overlaps to make the region complete
    isless(sum(overlapsAllowedHere), sum(missingPcHits[1:end - 1])) && return false, currentRegion
    # region is possible if it is an aggregate (i.e., no missing pcs)
    any(missingPcHits[1:end - 1]) || return true, currentRegion
    # region is not possible if not all missing pcs found in inPushPcs
    found, currentRegion = findMissingRegionPcs(currentRegion, aggregate, missingPcHits, overlapsAllowedHere) #  can ensure only overlaps corresponding to locations of missing are used?
    found || return false, currentRegion
    # region is not possible if not all duplicate pcs found in outPushPcs
    found, currentRegion = findDuplicateRegionPcs(currentRegion, overlapsAllowedHere) # find duplicate pcs COMBINE below into one or combine all with missingPcs
    found || return false, currentRegion
    return true, currentRegion
end

#= Find missing pcs from a region =#
"""
    findMissingRegionPcs()

Find the pcs not found in the region formed by a composition and determine whether or not these pcs are found in inPushPcs.  
"""
function findMissingRegionPcs(currentRegion, aggregate, missingPcHits, overlapsAllowedHere)
    currentRegion.missingPcs = aggregate[missingPcHits[1:end - 1]] # without 13's or duplicates
    # ensure that all missing pcs from the region are in inPushPcs
    for pc in currentRegion.missingPcs
        pc in currentRegion.inPushPcs[overlapsAllowedHere] || return false, currentRegion # if not in then return
    end
    return true, currentRegion
end

#= Find duplicate pcs from a region =#
"""
    findDuplicateRegionPcs()

Find the pcs found more than once in the region formed by a composition and determine whether or not these pcs are found in outPushPcs.    
"""
function findDuplicateRegionPcs(currentRegion, overlapsAllowedHere)
    duplicatePcs = zeros(Int8, length(currentRegion.missingPcs))
    duplicatePcHits = falses(length(currentRegion.pcs) + 1)
    # find duplciate pcs
    currentRegion.duplicatePcs = searchRegionPcsForDuplicatePcs(currentRegion, duplicatePcs, duplicatePcHits)
    # ensure that all duplicate pcs in the region are in outPushPcs
    testOutPushPcs = currentRegion.outPushPcs[overlapsAllowedHere]
    for pc in currentRegion.duplicatePcs
        dupFound = false
        for j::Int8 in 1:length(testOutPushPcs)
            if isequal(pc, testOutPushPcs[j])
                testOutPushPcs[j] = zero(eltype(testOutPushPcs))
                dupFound = true
                break
            end
        end
        dupFound || return false, currentRegion # if no dup found then return false
    end
    return true, currentRegion
end

"""
    searchRegionPcsForDuplicatePcs()

Determine the pcs found more than once in the region formed by a composition. Note that all occurences of 13 are duplicates. 
"""
function searchRegionPcsForDuplicatePcs(currentRegion, duplicatePcs, duplicatePcHits)
    k = one(Int8)
    for pc in currentRegion.pcs
        if duplicatePcHits[pc] || isless(12, pc) # 13
            duplicatePcs[k] = pc
            k += one(eltype(k))
        end
        duplicatePcHits[pc] = true
    end
    return duplicatePcs
end


#= Find all sets of overlaps allowed for a composition and the region it forms for a given position =#
"""
    findCandidateOverlapsForRegion()

Determine which (if any) possible overlaps for a composition successfully form a region containing an aggregate (i.e., 12 distinct pitch classes).    
"""
function findCandidateOverlapsForRegion(currentRegion, overlapsAllowedHere, allOverlaps, overlapIndicesByRow, rowsOverlapsByNum)
    # all sets of overlaps allowed for a composition at a given position and contain at least as many overlaps as the number of missing pcs 
    possibleOverlaps = findPossibleOverlapsForRegion(currentRegion, overlapsAllowedHere, allOverlaps, overlapIndicesByRow, rowsOverlapsByNum)
    indices = falses(size(possibleOverlaps, 1))
    for i in 1:size(possibleOverlaps, 1)    # search amtrix again instead?
        if isequal(sort(append!(getindex(currentRegion.inPushPcs, possibleOverlaps[i, :]), currentRegion.duplicatePcs)), sort(append!(getindex(currentRegion.outPushPcs, possibleOverlaps[i, :]), currentRegion.missingPcs)))
            setindex!(indices, true, i)
        end
    end
    return getindex(possibleOverlaps, indices, :) # error if only one overlap BitArray{1}?
end

"""
    findPossibleOverlapsForRegion()

Possible overlaps occur in matrix rows where (1) a composition contains a non-zero part, (2) each corresponding 
row position is greater than 1, and (3) there are at least as many as there are missingPcs. For example,
num. of missingPcs <= num. of possible overlaps <= num. of overlapsAllowedHere. This function reduces the 
number of overlaps that need to be tested.
"""
function findPossibleOverlapsForRegion(currentRegion, overlapsAllowedHere, allOverlaps, overlapIndicesByRow, rowsOverlapsByNum)   
    indices = trues(size(allOverlaps, 1))
    indices[overlapIndicesByRow[:, .!overlapsAllowedHere]] .= false # select only overlaps that share allowed rows
    if isless(1, length(currentRegion.missingPcs)) # select rows of overlaps with at least as many overlaps as the number of missing pcs
        indices[1:rowsOverlapsByNum[length(currentRegion.missingPcs)]] .= false
    end
    return getindex(allOverlaps, indices, :)
end  

#= Heuristics for guiding the search process by ranking candidate compostitions for a given position =#
"""
    selectHeuristicsAndCompute()

Selects and runs all heuristics chosen by the user.   
"""
function selectHeuristicsAndCompute(chosenHeuristics, position, idealPos, idealRegion, composition, candidateOverlaps, arraySpecs, prob)
    scoreSet = zeros(Float16, length(chosenHeuristics.heuristics))
    for (idx, choice) in enumerate(chosenHeuristics.heuristics)
        if isequal(choice, 0) # no heuristic
            return scoreSet, candidateOverlaps
        end
        if isequal(choice, 1)
            scoreSet[idx], candidateOverlaps = judgeQualityOfCandidate(idealRegion, candidateOverlaps, composition, position, idealPos)
        end
        # add additional heuristics here (would need to change user selection function as well)
    end
    return scoreSet, candidateOverlaps
end

"""
    judgeQualityOfCandidate()

An improved heuristic based on BemmanMeredithJNMR2016 where how far off a candidate is from an "ideal position" is
used to adjust the next "ideal position".
"""
function judgeQualityOfCandidate(idealRegion, candidateOverlaps, composition, position, idealPos)
    scores = zeros(Float16, size(candidateOverlaps,1))
    for i=1:length(scores)
        scores[i] = sum(abs, position + composition[1:end-1] - candidateOverlaps[i, :] .- 1 - idealPos)
        #scores[i] = abs(sum(abs, position + composition[1:end-1] - candidateOverlaps[i, :] .- 1 - idealPos) - idealRegion) # -1 necessary?
        #scores[i] = sum(position + composition[1:end-1] - candidateOverlaps[i, :] .- 1 - idealPos)
    end
    locs = sortperm(scores)
    candidateOverlaps = candidateOverlaps[locs,:]
    return scores[locs[1]], candidateOverlaps
end


#= Sort a list of candidates based on the ranked assigned to each by the heuristic(s) =#
"""
    sortCandidatesByHeuristicScore()

Order candidate compositions from best-to-worst according to the score assigned by the heuristic(s). 
Note that more than one heuristic can be used, in which case multiple scores are handled for a single candidate.
"""
function sortCandidatesByHeuristicScore(candidateScores)
    candidateScores = transpose(hcat(candidateScores...)) # convert to matrix
    # normalize scores
    for i in 1:size(candidateScores,2)
        candidateScores[:,i] = (candidateScores[:,i] .- minimum(candidateScores[:,i]))/(maximum(candidateScores[:,i]) - minimum(candidateScores[:,]))
    end
    return sortperm(vec(sum(candidateScores, dims=2))) # single vector of summed normalized heuristic scores
end



#= Advancing the search when a candidate has been selected for a given position =#
"""
    success()

Proceeds after selecting the current candidate.
"""
function success(position, cnt, chosenCandidate, overlapIdx)
    position += chosenCandidate.composition - chosenCandidate.overlaps[overlapIdx, :]
    cnt += one(eltype(cnt))
    return cnt, position, chosenCandidate.partition
end

#= Backtracking the search when no more candidates can be selected for a given position =#
"""
    backtrack()

Backtracks after deselecting the current candidate.
"""
function backtrack(position, cnt, chosenCandidate, overlapIdx)
    cnt -= one(eltype(cnt))
    position -= chosenCandidate.composition - chosenCandidate.overlaps[overlapIdx, :]
    return cnt, position, Array{Candidate,1}(), zero(Int8), zero(Int32), zero(Int16)
end

function backtrack(cnt)
    cnt -= one(eltype(cnt))
    return cnt, Array{Candidate,1}(), zero(Int8), zero(Int32), zero(Int16)
end

function backtrack(cnt, position, allCandidateLists, usedPartitions, selected)
    allCandidateLists[cnt] = Array{Candidate,1}()
    usedPartitions[cnt] = zero(Int8)
    selected.candidates[cnt] = zero(Int32) # use only candidates to control (allCandidatesLists is too big)?
    selected.overlaps[cnt] = zero(Int16)
    cnt -= one(eltype(cnt))
    position -= allCandidateLists[cnt][selected.candidates[cnt]].composition - allCandidateLists[cnt][selected.candidates[cnt]].overlaps[selected.overlaps[cnt], :]
    return cnt, position, allCandidateLists, usedPartitions, selected
end

#= Prints the furthest progress reached in the search =#
"""
    printProgress()

Prints to console the current status of the search.
"""
function printProgress(selected, farthestReached, cnt)
    # >= prints all partial solutions at least as good as the farthest so far
    # > print only the farthest so far
    if cnt > farthestReached
        farthestReached = cnt
        #println(selected.candidates)
    end
    return farthestReached
end

#= Prints the first solution found =#
"""
    printSolution()

Prints to console the results of the search.
"""
function printSearchResults(selected, allCandidateLists, cnt, inputMatrix, numBacktracks)
    if selected.candidates[1] == 0
        println("No solution possible.")
    end
    #println("All successful candidates with all overlaps:")
    #println("Num.", "\t", "Part.", "\t", "Composition", "\t", "\t", "\t", "Overlaps")
    #for i = 1:cnt - 1
    #    println(i, "\t", allCandidateLists[i][cvec[i]].partition, "\t", allCandidateLists[i][cvec[i]].composition, "\t", "\t", "\t", allCandidateLists[i][cvec[i]].overlaps)
    #end
    println("\n", "Number of backtracks required: ", numBacktracks)
    println("\n", "All successful candidates and their chosen overlaps in this solution:", "\n")
    println("Num.", "\t", "Part.", "\t", "Composition", "\t", "\t", "\t", "Overlaps")
    for i = 1:cnt - 1
        println(i, "\t", allCandidateLists[i][selected.candidates[i]].partition, "\t", allCandidateLists[i][selected.candidates[i]].composition, "\t", "\t", allCandidateLists[i][selected.candidates[i]].overlaps[selected.overlaps[i],:])
    end
    println("\n", "All regions formed by the successful candidates in this solution:")
    position = ones(Int64, length(allCandidateLists[1][selected.candidates[1]].composition))
    for i = 1:cnt - 1
        println("\n", "Region ", i)
        candidateComp = allCandidateLists[i][selected.candidates[i]].composition
        candidateOverlaps = allCandidateLists[i][selected.candidates[i]].overlaps[selected.overlaps[i],:]
        for j=1:length(position)
            startIdx = position[j] - candidateOverlaps[j]
            endIdx = startIdx + candidateComp[j] - 1
            if candidateComp[j] > 0
                println(j, " ", inputMatrix[j, startIdx:endIdx])
            else
                println(j)
            end
        end
        position = position + candidateComp - candidateOverlaps
    end
    return nothing
end

#= End of file =#
