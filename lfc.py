# Solve Manning's equation in order to test for the low-flow channel thing.

from time import time
import sys
DEBUG = False

lfcDepth = 1
lfcWidth = 10
cWidth = 30 # Incl. LFC
cDepth = 5 # above LFC
slope = 0.001
c = 1.49
n = 0.01

# Estimate based on usage (running on hydro lab computer)
# 120k iterations in 10 seconds
TYPICAL_RATE = 10.0 / 100000

def hasJump(data, method = "local", factor = 3, localRange = 5, relativeFactor = 0.05):
    if DEBUG:
        print(method)
    # Determine if there is a sudden jump in a set of points over equal intervals
    # Assume that direction of change is monotonic
    # Factor: step must exceed factor * reference step to be counted as jump.
    # localRange: number of steps to compare for "local" method
    # relativeFactor: maximum value for step / previous value
    if method == "average":
        # Average: compare step at each pair of elements to overall average step.
        # Note that this may cause issues if the data is very non-linear,
        # as the average step may be very high or very low compared to a local area
        if len(data) < 2:
            return False # If there is only 1 element, or none, there can't be a jump
        else:
            avgStep = abs(float(data[-1] - data[0]) / (len(data) - 1)) # Needs to be float so we don't end up with 0 if the step is small
            foundJump = False
            ix = 0
            while not (foundJump or ix >= len(data) - 1):
                delta = data[ix + 1] - data[ix]
                if abs(delta) > factor * avgStep:
                    foundJump = True
                    return True
                ix = ix + 1
            return foundJump
    if method == "local":
        # Local: compare each step to the preceding and next steps
        # This might not find a jump if it is split over multiple steps.
        if len(data) < 2:
            return False # 1 element -> no jump
        if len(data) < (localRange + 1):
            raise ValueError("Local method needs at least %d steps (%d elements) to work" % (localRange, localRange + 1))
        ix = 0
        foundJump = False
        while not (foundJump or ix > len(data) - (localRange + 1)):
            testRange = data[ix:(ix+localRange+1)]
            # Use average only on a local range
            foundJump = hasJump(testRange, method = "average", factor = factor)
            if foundJump:
                return True
            ix = ix + 1
        return foundJump
    if method == "relative":
        # Relative: compare each step to the depth at the previous point.
        # This won't work on very steep curves, but it may be more accurate
        # to the general intent of the calibrator.
        # If a relative jump is found, it will also test local.
        if len(data) < 2:
            return False
        ix = 0
        foundJump = False
        while not (foundJump or ix >= len(data) - 1):
            foundJump = (abs(float(data[ix + 1] - data[ix]) / data[ix]) > relativeFactor) and hasJump(data, "local", factor, localRange)
            if DEBUG:
                print("%f, %f: %f" % (data[ix+1], data[ix], abs((data[ix+1] - data[ix]) / data[ix])))
            if foundJump:
                return True
            ix = ix + 1
        return foundJump
    # Return from each "if" statement -- the code should only get here if
    # no method is matched
    raise ValueError("Error: invalid method specified.  Method must be 'average', 'local', or 'relative'.")

def area(depth, lfcDepth, lfcWidth, cWidth, z):
    # z = 0 for rectangular
    topWidth = lfcWidth + 2 * lfcDepth * z if depth >= lfcDepth else lfcWidth + 2 * depth * z
    if depth < lfcDepth:
        return depth * (lfcWidth + topWidth) / 2
    else:
        return (lfcDepth * (lfcWidth + topWidth) / 2) + (depth - lfcDepth) * cWidth
        
def wp(depth, lfcDepth, lfcWidth, cWidth, z):
    # z = 0 for rectangular
    if depth < lfcDepth:
        return lfcWidth + 2 * depth * (z**2 + 1)**(1/2)
    else:
        return lfcWidth + 2 * lfcDepth * (z**2 + 1)**(1/2) + (cWidth - lfcWidth) + 2 * (depth - lfcDepth)
        
def radius(depth, lfcDepth, lfcWidth, cWidth, z):
    return area(depth, lfcDepth, lfcWidth, cWidth, z) / wp(depth, lfcDepth, lfcWidth, cWidth, z)
    
def manningQ(depth, c, lfcn, mcn, lfcDepth, lfcWidth, cWidth, slope, z):
    # mcn: main channel n
    # Calculate WP-weighted n
    if depth == 0:
        return 0
    totalWp = wp(depth, lfcDepth, lfcWidth, cWidth, z)
    # This produces the maximum lfcWp, but if the depth < lfcDepth then the result will still just be lfcn
    lfcWp = wp(lfcDepth, lfcDepth, lfcWidth, cWidth, z)
    mcWp = totalWp - lfcWp if totalWp > lfcWp else 0
    compositeN = (lfcn * lfcWp + mcn * mcWp) / totalWp
    if DEBUG and compositeN == 0:
        print("CompositeN: %f totalWp: %f lfcWp: %f mcWp: %f lfcn: %f mcn: %f" % (compositeN, totalWp, lfcWp, mcWp, lfcn, mcn))
    return c / compositeN * area(depth, lfcDepth, lfcWidth, cWidth, z) * radius(depth, lfcDepth, lfcWidth, cWidth, z)**(2/3) * slope ** (1/2)
    
def solveDepth(q, lfcn, mcn, lfcDepth, lfcWidth, cWidth, slope, z, c = 1.49, tolerance = 0.01, lutolerance = 0.0000000001):
    # Sort of a binary search thing
    # Initial bounds are arbitrary, but lower needs to be > 0
    # lutolerance is set in order to avoid infinite loop when lower and upper
    # become very close together but can't solve (why does that happen?)
    lower = 0.01
    upper = 100.0
    depth = 0
    # First, find what the bounds should actually be
    while manningQ(lower, c, lfcn, mcn, lfcDepth, lfcWidth, cWidth, slope, z) > q:
        lower = lower / 2
    while manningQ(upper, c, lfcn, mcn, lfcDepth, lfcWidth, cWidth, slope, z) < q:
        upper = upper * 2
    
    found = False
    while not found:
        depth = (lower + upper) / 2
        mq = manningQ(depth, c, lfcn, mcn, lfcDepth, lfcWidth, cWidth, slope, z)
        if abs(mq - q) <= tolerance:
            found = True
            break
        if mq > q:
            upper = depth
        if mq < q:
            lower = depth
        if lower + lutolerance >= upper: # Couldn't find depth at all
            # Assume that this is the closest we can get?
            found = True
    if found:
        return depth
    else:
        return False
            
            
def testPerms(qs, lfcns, mcnfs, lfcds, lfcws, cws, slopes, zs, c = 1.49, widthFactor = 3, jumpMethod = "local", jumpFactor = 3, relativeFactor = 0.05, log = True, increments = 20):
    # Test all permutations within the given ranges for whether they have a jump in the depths
    # If log, then notify the user every 1/increments of the total number of iterations
    # Also notify the user if a no-jump combination has been found.
    # widthFactor: minimum value for cw / lfcw in order to test combination
    # mcnfs: Factors for mcn / lfcn, instead of testing unrelated sets of ns
    depthsets = []
    perms = len(ns) * len(lfcds) * len(lfcws) * len(cws) * len(slopes) * len(zs) * len(mcnfs)
    iters = perms * len(qs)
    # Counting is easier with permutations, but timing should be done with iterations since this
    # accounts for the number of flow rates to be tested
    incrSize = perms // increments + 1
    if log:
        print("Max total permutations: %d (iterations: %d)" % (perms, iters))
        t = TYPICAL_RATE * iters // 60
        units = "minutes"
        if t == 0:
            t = TYPICAL_RATE * iters
            units = "seconds"
        print("**VERY** rough time estimate (typically an overestimate): %d %s" % (t, units))
    start = time()
    for n in ns:
        for lfcd in lfcds:
            for lfcw in lfcws:
                for cw in cws:
                    if cw >= lfcw * widthFactor: # It is absurd for the main channel to be narrower than the LFC, and only useful if it is significantly wider
                        for slope in slopes:
                            if slope > 0: # 0 slope simply won't flow
                                for z in zs:
                                    if (lfcw + 2 * z * lfcd < (cw / widthFactor)): # Can't have LFC top width > cw
                                        for mcnf in mcnfs:
                                            depths = [solveDepth(q, n, n * mcnf, lfcd, lfcw, cw, slope, z) for q in qs]
                                            # Sometimes q is in the middle of a jump; remove that
                                            depths = [h for h in depths if h]
                                            jump = hasJump(depths, jumpMethod, jumpFactor, relativeFactor = relativeFactor)
                                            depthsets.append({
                                                "jump": jump,
                                                "n": n,
                                                "mcnf": mcnf,
                                                "lfcd": lfcd,
                                                "lfcw": lfcw,
                                                "cw": cw,
                                                "slope": slope,
                                                "z": z,
                                                "depths": depths
                                                })
                                            if log:
                                                """if not jump:
                                                    print("No-jump combination found with n = %f, lfcd = %f, lfcw = %f, cw = %f (max depth = %f)" % (n, lfcd, lfcw, cw, depths[-1]))"""
                                                ld = len(depthsets)
                                                if ld % incrSize == 0:
                                                    pctDone = int(ld * len(qs) * 100 / iters)
                                                    # Base length is 54, progress bar length is 27
                                                    string = "Progress: completed %d iterations (%d%%) in %d seconds" % (ld * len(qs), pctDone, time() - start)
                                                    progressBar = (" " * (73 - len(string) if len(string) < 72 else 1)) + "[%s]" % ("#" * (pctDone // 5) + " " * (20 - pctDone // 5))
                                                    print(string + progressBar)
    return depthsets
    
def depthsetString(depthsets, prt = True, table = True):
    # Format depthsets in a structured way so patterns are visible
    # Extra check to make sure the no-jump isn't just because it never leaves the LFC or is never in it
    PRINTWIDTH = 15
    depthsets = [p for p in depthsets if not p["jump"]
        and (p["depths"][-1] > 1.1*p["lfcd"])
        and (p["depths"][0] < p["lfcd"])]
    uniqueNs = list(set([p["n"] for p in depthsets]))
    out = []
    outStr = ""
    if not table: # Print in an expanded format - more human-readable
        out.append("Parameter sets with no jump:")
        for n in uniqueNs:
            out.append("n: %f" % n)
            dsn = [p for p in depthsets if p["n"] == n]
            uniqueLfcds = list(set([p["lfcd"] for p in dsn]))
            for lfcd in uniqueLfcds:
                out.append("\tLFCD: %f" % lfcd)
                dslfcd = [p for p in dsn if p["lfcd"] == lfcd]
                uniqueLfcws = list(set([p["lfcw"] for p in dslfcd]))
                for lfcw in uniqueLfcws:
                    out.append("\t\LFCW: %f" % lfcw)
                    dslfcw = [p for p in dsn if p["lfcw"] == lfcw]
                    uniqueCws = list(set([p["cw"] for p in dslfcw]))
                    for cw in uniqueCws:
                        out.append("\tCW: %f" % cw)
                        dscw = [p for p in dslfcw if p["cw"] == cw]
                        uniqueSlopes = list(set([p["slope"] for p in dscw]))
                        for slope in uniqueSlopes:
                            out.append("\tSlope: %f" % slope)
                            dsslope = [p for p in dscw if p["slope" == slope]]
                            uniqueZs = list(set([p["z"] for p in dsslope]))
                            for z in uniqueZs:
                                out.append("\tz: %f" % z)
        outStr = "\n".join(out)
    if table: # Make a table format
        order = ["n", "MC n Factor", "LFC Depth (ft)", "LFC Width (ft)", "MC Width (ft)", "Slope", "z", "Min Depth (ft)", "Max Depth (ft)"]
        paramsets = [[p["n"], p["mcnf"], p["lfcd"], p["lfcw"], p["cw"], p["slope"], p["z"], p["depths"][0], p["depths"][-1]] for p in depthsets]
        paramsets.sort()
        out = [order] + paramsets
        outStr = "\n".join(["\t".join([" " + i[:15] if len(i) > 15 else (" " * (16 - len(i)) + i) for i in [str(i) for i in l]]) for l in out])
    if prt:
        print(outStr)
    return out
    
def findOptimal(qs, lfcns, mcnfs, lfcds, lfcws, cws, slopes, zs, jumpMethod = "local", jumpFactor = 3, log = True, upperBoundJ = 5, upperBoundR = 2, tolerance = 0.01, iterLimit = 20, printSets = True):
    # Use a "binary search" method to find the optimal combination of parameters
    # Tolerance is the precision of bounds finding
    # Upper bound: maximum jump or relative factor to assume that no reasonable solution will be found
    upper = float(upperBoundR if jumpMethod == "relative" else upperBoundJ)
    lower = float(0 if jumpMethod == "relative" else 1) # The largest step will always be >= the average step
    iterResults = []
    def iteration(factor):
        # Find how many no-jump results there are for a given factor
        jf = jumpFactor if jumpMethod == "relative" else factor
        rf = factor
        depthsets = testPerms(qs, lfcns, mcnfs, lfcds, lfcws, cws, slopes, zs, jumpMethod = jumpMethod, jumpFactor = jf, relativeFactor = rf, log = log)
        depthsets = [p for p in depthsets if not p["jump"]
            and (p["depths"][-1] > 1.1*p["lfcd"])
            and (p["depths"][0] < p["lfcd"])]
        if len(depthsets) > 0:
            # Only set iter results if there are results
            iterResults = depthsets
        if printSets:
            depthsetString(depthsets)
        return len(depthsets)
    # First, check that the upper bound isn't exceeded
    if (iteration(upper) == 0):
        if log:
            print("No no-jump solution found with method %s" % method)
        return False
    else:
        # Keep iterating until a single (within tolerance) point is found where solutions = 1
        count = 0
        while upper > (lower + tolerance) and count < iterLimit: # Don't go on running forever
            factor = (upper + lower) / 2
            if log:
                print("Testing with factor %f" % factor)
            solutions = iteration(factor)
            if log:
                print("Found %d solutions" % solutions)
            if solutions > 0: # Current factor is an upper bound
                upper = factor
            if solutions == 0: # Current factor is a lower bound
                lower = factor
            count = count + 1
        if log:
            print("Found upper bound %f (lower bound %f, iterations %d)" % (upper, lower, count))
        return (upper, iterResults)
    
def writeCsv(depthsetTable, path):
    print("Writing CSV to %s" % path)
    with open(path, "w") as of:
        of.write("\n".join([",".join([str(f) for f in i]) for i in depthsetTable]))
        
def ratingCurve(qs, lfcn, mcnf, lfcDepth, lfcWidth, cWidth, slope, z, graphWidth = 25):
    # Print out a rating curve graph to the console
    mcn = lfcn * mcnf
    depths = [(q, solveDepth(q, lfcn, mcn, lfcDepth, lfcWidth, cWidth, slope, z)) for q in qs]
    maxDepth = depths[-1][1]
    increment = float(maxDepth) / graphWidth # Value per graph character
    def numString(q, h):
        sq = str(q)
        sh = str(h)
        return " " * (10 - len(sq) if len(sq) <= 10 else 0) + \
            (sq[:10] if len(sq) > 10 else sq) + " : " + \
            (sh[:10] if len(sh) > 10 else sh)
    for (q, h) in depths:
        sout = "%s %s" % (numString(q, h), "X" * int(h / increment))
        print(sout)
    
helpStr = """Low Flow Channel Permutations:
This program tests if there are any permutations of parameters
which do not result in an LFC jump.  Low-flow and main channels
are assumed rectangular for now.

Note that progress will be underestimated, as it is based on the maximum
total permutations, whereas fewer than that actually happen due to filtering
out absurd combinations.

Parameters:
Format: -identifier<param1-p2-...-pn>, e.g. -q1-5-10
All <min-max-step> parameters can be replaced by <value> in order to only test one value.

    -q<min-max-step>:       Minimum, maximum, and step flow rates (cfs).  Default 1, 100, 10.
    -n<min-max-step>:       Same thing, but Manning's n.  Default 0.01, 0.1, 0.01.
    -m<min-max-step>:       Main channel n factors, i.e. main channel n / lfc n.  Default 1, 1, 1.
    -d<min-max-step>:       LFC depths (ft).  Default 0.5, 10, 1.
    -w<min-max-step>:       LFC widths (ft).  Default 1, 50, 5.
    -c<min-max-step>:       Channel widths (ft).  Default 10, 100, 10.
    -j<factor>      :       Factor for jump detection.  Default 3.
    -f[path]        :       Write output to a CSV file at <path>.  Path is optional, but the default is specific
                                to the developer's use case, so other users should always set the path.  The default
                                is "Z:\\adit\\Desktop\\LARFlows\\Data Processing\\Calibration\\lfc.csv".
    -a                      Use "average" jump detection method instead of "local".
    -z<min-max-step>:       z for trapezoidal LFC.  Default 0, 1, 1, i.e. 0 - rectangle.
    -s<min-max-step>:       Slope.  Default 0.001, 0.001, 0.001.
    -r[factor]      :       Use "relative" jump detection method with, optionally, relative factor = factor.
                                The default factor is 0.05.
    -p              :       Print out rating curve.  Will use the first value in each range except for q.
    -o              :       Search for the optimal solution under the given parameters.  Will run many iterations
                                and likely take quite some time.  In this case -j specifies jump factor for use
                                with relative method and -r is ignored.
    -h              :       Show this help message.
"""

def floatRange(start, end, step):
    # Because built-in range doesn't support floats
    while start <= end:
        yield start
        start += step

if __name__ == "__main__":
    qs = list(floatRange(1, 1000, 10))
    ns = list(floatRange(0.01, 0.1, 0.01))
    mcnfs = [1]
    lfcds = list(floatRange(0.5, 10, 1))
    lfcws = list(floatRange(1, 50, 5))
    cws = list(floatRange(10, 100, 10))
    method = "local"
    rf = 0.05
    jf = 3
    zs = [0]
    slopes = [0.001]
    write = False
    path = "Z:\\adit\\Desktop\\LARFlows\\Data Processing\\Calibration\\lfc.csv"
    args = sys.argv[1:]
    help = False
    rc = False
    optimal = False
    for arg in args:
        if arg == "-h" or arg == "--help" or arg == "help":
            help = True
            break
        if len(arg) >= 2:
            specifier = arg[1]
            vals = arg[2:].split("-")   # List of argument values
            rg = []                     # Range from arguments
            num = 0                     # Numerical argument
            st = arg[2:]                # String argument
            if specifier in "qndwcszm": # Range arguments
                if len(vals) == 3:
                    vals = [float(a) for a in vals]
                    rg = list(floatRange(vals[0], vals[1], vals[2]))
                if len(vals) == 1 and vals[0] != "":
                    rg = [float(vals[0])]
                if not len(vals) in [1, 3]:
                    print("Error: argument %s requires either 1 or 3 values." % specifier)
                    help = True
                    break
            if specifier in "j":        # Number arguments
                num = float(st)
            if specifier == "q":
                qs = rg
            if specifier == "n":
                ns = rg
            if specifier == "d":
                lfcds = rg
            if specifier == "w":
                lfcws = rg
            if specifier == "c":
                cws = rg
            if specifier == "j":
                jf = num
            if specifier == "z":
                zs = rg
            if specifier == "s":
                slopes = rg
            if specifier == "m":
                mcnfs = rg
            if specifier == "o":
                optimal = True
            if specifier == "r":
                method = "relative"
                if len(vals) == 1 and vals[0] != "":
                    rf = float(vals[0])
            if specifier == "f":
                write = True
                if st != "":
                    path = st
            if specifier == "a":
                method = "average"
            if specifier == "p":
                rc = True
            if not specifier in "qndwcjfaszmrpo":
                print("Error: nonexistent parameter specified.")
                help = True
                break
        if len(arg) <= 2 and not (arg in ["-h", "-f", "-a", "-r", "-p", "-o"]):
            print("Error: nonexistent parameter specified or parameter values not provided.")
            help = True
            break
    if help:
        print(helpStr)
    elif rc:
        ratingCurve(qs, ns[0], mcnfs[0], lfcds[0], lfcws[0], cws[0], slopes[0], zs[0])
    elif optimal:
        (jOpt, jIter) = findOptimal(qs, ns, mcnfs, lfcds, lfcws, cws, slopes, zs, jumpMethod = "local")
        (rOpt, rIter) = findOptimal(qs, ns, mcnfs, lfcds, lfcws, cws, slopes, zs, jumpMethod = "relative", jumpFactor = jf)
        print("Optimal j: %f; optimal r: %f" % (jOpt, rOpt))
        print("Optimal j parameters:")
        depthsetString(jIter)
        print("Optimal r parameters:")
        depthsetString(rIter)
    else:
        result = depthsetString(testPerms(qs, ns, mcnfs, lfcds, lfcws, cws, slopes, zs, jumpFactor = jf, jumpMethod = method, relativeFactor = rf))
        if write:
            writeCsv(result, path)