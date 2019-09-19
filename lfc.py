# Solve Manning's equation in order to test for the low-flow channel thing.

from time import time
import sys

lfcDepth = 1
lfcWidth = 10
cWidth = 30 # Incl. LFC
cDepth = 5 # above LFC
slope = 0.001
c = 1.49
n = 0.01

# Estimate based on usage (running on hydro lab computer)
# 17000 iterations in 514 seconds
TYPICAL_RATE = 500.0 / 17000

def hasJump(data, method = "average", factor = 3, localRange = 5):
    # Determine if there is a sudden jump in a set of points over equal intervals
    # Assume that direction of change is monotonic
    # Factor: step must exceed factor * reference step to be counted as jump.
    # localRange: number of steps to compare for "local" method
    if method == "average":
        # Average: compare step at each pair of elements to overall average step.
        # Note that this may cause issues if the data is very non-linear,
        # as the average step may be very high or very low compared to a local area
        if len(data) < 2:
            return False # If there is only 1 element, or none, there can't be a jump
        else:
            avgStep = float(data[-1] - data[0]) / (len(data) - 1) # Needs to be float so we don't end up with 0 if the step is small
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
    # Return from each "if" statement -- the code should only get here if
    # no method is matched
    raise ValueError

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
    totalWp = wp(depth, lfcDepth, lfcWidth, cWidth, z)
    lfcWp = wp(lfcDepth, lfcDepth, lfcWidth, cWidth, z)
    mcWp = totalWp - lfcWp
    compositeN = (lfcn * lfcWp + mcn * mcWp) / totalWp
    return c / compositeN * area(depth, lfcDepth, lfcWidth, cWidth, z) * radius(depth, lfcDepth, lfcWidth, cWidth, z)**(2/3) * slope ** (1/2)
    
def solveDepth(q, lfcn, mcn, lfcDepth, lfcWidth, cWidth, slope, z, c = 1.49, tolerance = 0.01):
    # Sort of a binary search thing
    # Initial bounds are arbitrary, but lower needs to be > 0
    lower = 0.01
    upper = 100
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
        if (mq < q + tolerance) and (mq > q - tolerance):
            found = True
            break
        if mq > q:
            upper = depth
        if mq < q:
            lower = depth
        if lower >= upper: # Couldn't find depth at all
            break
    if found:
        return depth
    else:
        return False
            
            
def testPerms(qs, lfcns, mcnfs, lfcds, lfcws, cws, slopes, zs, c = 1.49, widthFactor = 3, jumpMethod = "local", jumpFactor = 3, log = True, increments = 100):
    # Test all permutations within the given ranges for whether they have a jump in the depths
    # If log, then notify the user every 1/increments of the total number of iterations
    # Also notify the user if a no-jump combination has been found.
    # widthFactor: minimum value for cw / lfcw in order to test combination
    # mcnfs: Factors for mcn / lfcn, instead of testing unrelated sets of ns
    depthsets = []
    iters = len(ns) * len(lfcds) * len(lfcws) * len(cws) * len(slopes) * len(zs) * len(mcnfs)
    if log:
        print("Max total permutations: %d" % iters)
        print("**VERY** rough time estimate: %d minutes" % (TYPICAL_RATE * iters / 60))
    start = time()
    for n in ns:
        for lfcd in lfcds:
            for lfcw in lfcws:
                for cw in cws:
                    if cw >= lfcw * widthFactor: # It is absurd for the main channel to be narrower than the LFC, and only useful if it is significantly wider
                        for slope in slopes:
                            if slope > 0: # 0 slope simply won't flow
                                for z in zs:
                                    for mcnf in mcnfs:
                                        depths = [solveDepth(q, n, n * mcnf, lfcd, lfcw, cw, slope, z) for q in qs]
                                        jump = hasJump(depths, jumpMethod, jumpFactor)
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
                                            if (ld * increments) % iters == 0:
                                                pctDone = int(ld * 100 / iters)
                                                print("Progress: completed %d iterations (%d%%) in %d seconds \t[%s]" % (ld, pctDone, time() - start, ("#" * (pctDone // 4) + " " * ((100 - pctDone) // 4))))
    return depthsets
    
def depthsetString(depthsets, prt = True, table = True):
    # Format depthsets in a structured way so patterns are visible
    # Extra check to make sure the no-jump isn't just because it never leaves the LFC or is never in it
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
        outStr = "\n".join(["\t\t".join([i[:15] if len(i) > 15 else i for i in [str(i) for i in l]]) for l in out])
    if prt:
        print(outStr)
    return out
    
def writeCsv(depthsetTable, path):
    print("Writing CSV to %s" % path)
    with open(path, "w") as of:
        of.write("\n".join([",".join([str(f) for f in i]) for i in depthsetTable]))

    
helpStr = """Low Flow Channel Permutations:
This program tests if there are any permutations of parameters
which do not result in an LFC jump.  Low-flow and main channels
are assumed rectangular for now.

Note that progress will be underestimated, as it is based on the maximum
total permutations, whereas fewer than that actually happen due to filtering
out absurd combinations.

Parameters:
Format: -identifier<param1-p2-...-pn>, e.g. -q1-5-10

    -q<min-max-step>:       Minimum, maximum, and step flow rates (cfs).  Default 1, 100, 10.
    -n<min-max-step>:       Same thing, but Manning's n.  Default 0.01, 0.1, 0.01.
    -m<min-max-step>:       Main channel n factors, i.e. main channel n / lfc n.  Default 1, 2, 1, i.e. 1.
    -d<min-max-step>:       LFC depths (ft).  Default 0.5, 10, 1.
    -w<min-max-step>:       LFC widths (ft).  Default 1, 50, 5.
    -c<min-max-step>:       Channel widths (ft).  Default 10, 100, 10.
    -j<factor>      :       Factor for jump detection.  Default 3.
    -f[path]        :       Write output to a CSV file at <path>.  Path is optional, but the default is specific
                                to the developer's use case, so other users should always set the path.
    -a                      Use "average" jump detection method instead of "local".
    -z<min-max-step>        z for trapezoidal LFC.  Default 0, 1, 1, i.e. 0 - rectangle.
    -s<min-max-step>        Slope.  Default 0.001, 0.002, 0.001, i.e. 0.001.
    -h              :       Show this help message.
"""

def floatRange(start, end, step):
    # Because built-in range doesn't support floats
    while start < end:
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
    jf = 3
    zs = [0]
    slopes = [0.001]
    write = False
    path = "Z:\\adit\\Desktop\\LARFlows\\Data Processing\\Calibration\\lfc.csv"
    args = sys.argv[1:]
    help = False
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
                if len(vals) != 3:
                    print("Error: parameters to -%s require 3 values." % specifier)
                    help = True
                    break
                vals = [float(a) for a in vals]
                rg = list(floatRange(vals[0], vals[1], vals[2]))
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
            if specifier == "f":
                write = True
                if st != "":
                    path = st
            if specifier == "a":
                method = "average"
            if not specifier in "qndwcjfaszm":
                print("Error: nonexistent parameter specified.")
                help = True
                break
        if len(arg) <= 2 and not (arg in ["-h", "-f", "-a"]):
            print("Error: nonexistent parameter specified or parameter values not provided.")
            help = True
            break
    if help:
        print(helpStr)
    if not help:
        result = depthsetString(testPerms(qs, ns, mcnfs, lfcds, lfcws, cws, slopes, zs, jumpFactor = jf))
        if write:
            writeCsv(result, path)