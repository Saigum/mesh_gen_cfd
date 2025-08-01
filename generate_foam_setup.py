import os
import textwrap
from PyFoam.Execution.BasicRunner import BasicRunner
import os


FOAM_BASHRC = "/usr/lib/openfoam/openfoam2412/etc/bashrc"

def run_with_basic_runner(cmd):
    """
    Run an OpenFOAM utility with PyFoam's BasicRunner,
    after sourcing the OpenFOAM-2412 environment.
    """
    # build a single shell command that sources foam and then runs the utility
    shell_cmd = f"source /usr/lib/openfoam/openfoam2412/etc/bashrc && {' '.join(cmd)}"
    print(shell_cmd)
    runner = BasicRunner(
        argv=["bash", "-lc", shell_cmd],
        silent=False  
    )
    runner.start()
    
def mesh_generation(
    # ─────────────────────────────────────────────────────────────────────────────
    # 1) blockMeshDict parameters
    convert_to_meters: float = 1.0,
    vertices: tuple = (
        # (xmin, ymin, zmin)
        ( 40.6494,   76.7822,   18.750),
        # (xmax, ymin, zmin)
        (161.694,    76.7822,   18.750),
        # (xmax, ymax, zmin)
        (161.694,   168.560,    18.750),
        # (xmin, ymax, zmin)
        ( 40.6494,  168.560,    18.750),
        # (xmin, ymin, zmax)
        ( 40.6494,   76.7822,  104.250),
        # (xmax, ymin, zmax)
        (161.694,    76.7822,  104.250),
        # (xmax, ymax, zmax)
        (161.694,   168.560,   104.250),
        # (xmin, ymax, zmax)
        ( 40.6494,  168.560,   104.250),
    ),
    block_size: tuple = (20, 20, 20),
    simple_grading: tuple = (1, 1, 1),

    # ─────────────────────────────────────────────────────────────────────────────
    # 2) createPatchDict parameters
    inlet_box: tuple = (
        (-1e5, -1e5, -1e5),
        ( 41.15,  1e5,  1e5),
    ),
    outlet_box: tuple = (
        (161.20, -1e5, -1e5),
        ( 1e5,   1e5,  1e5),
    ),

    # ─────────────────────────────────────────────────────────────────────────────
    # 3) controlDict parameters
    application: str = "snappyHexMesh",
    startTime: float = 0.0,
    endTime: float = 1.0,       # dummy, only a mesh run
    deltaT: float = 1.0,
    writeInterval: int = 1,

    # ─────────────────────────────────────────────────────────────────────────────
    # 4) fvSchemes – leave defaults as shown
    ddtScheme: str = "Euler",
    gradScheme: str = "Gauss linear",
    divScheme: str = "none",
    laplacianScheme: str = "Gauss linear corrected",
    interpolationScheme: str = "linear",
    snGradScheme: str = "corrected",

    # ─────────────────────────────────────────────────────────────────────────────
    # 5) meshQualityDict – defaults as given
    maxNonOrtho: int         = 65,
    maxBoundarySkewness: int = 20,
    maxInternalSkewness: int = 4,
    maxConcave: int          = 80,
    minFlatness: float       = 0.5,
    minVol: float            = 1e-13,
    minTetQuality: float     = 1e-30,
    minArea: float           = -1.0,
    minTwist: float          = 0.02,
    minDeterminant: float    = 0.001,
    minFaceWeight: float     = 0.02,
    minVolRatio: float       = 0.01,
    minTriangleTwist: float  = -1.0,
    nSmoothScale: int        = 4,
    errorReduction: float    = 0.75,

    # ─────────────────────────────────────────────────────────────────────────────
    # 6) snappyHexMeshDict parameters
    castellatedMesh: bool    = True,
    snap: bool               = True,
    addLayers: bool          = False,
    stl_filename: str        = "output_surface.stl",  # must live in constant/triSurface/
    refinement_levels: tuple = (1, 2),
    maxLocalCells: int       = 200_000,
    maxGlobalCells: int      = 500_000,
    nCellsBetweenLevels: int = 2,
    locationInMesh: tuple    = (101.17, 122.67,  61.50),
    resolveFeatureAngle: int = 30,
    mergeTolerance: float    = 1e-6,

    # ─────────────────────────────────────────────────────────────────────────────
    # 7) surfaceFeatureExtractDict parameters
    input_surface: str = "output_surface.stl",
    triContinuity: int = 1,
    writeObj: bool    = True,
    writeSTL: bool    = False,
    writeThree: bool  = False,
    writeASCII: bool  = True,
    prefix: str = "vessel",
):
    """
    Create all necessary system/ and constant/ dictionaries for a standard
    snappyHexMesh workflow (assuming you already have your STL in constant/triSurface/).
    """

    
    ## this is for the mesh generation step, so i'm making the openfoam setup here in a directory called mesh
    
    # 1) Make sure directories exist
    os.makedirs(f"{prefix}/mesh/system", exist_ok=True)
    os.makedirs(f"{prefix}/mesh/constant", exist_ok=True)
    os.makedirs(f"{prefix}/mesh/constant/polySurface", exist_ok=True)
    os.makedirs(f"{prefix}/mesh/constant/triSurface", exist_ok=True)
    os.makedirs(f"{prefix}/mesh/0", exist_ok=True)

    # ─────────────────────────────────────────────────────────────────────────────
    # 2) Write blockMeshDict
    bm_content = f"""
    /*--------------------------------*- C++ -*----------------------------------*\\
    | =========                 |                                                 |
    | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\\\    /   O peration     | Version:  2412                                  |
    |   \\\\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\\\/     M anipulation  |                                                 |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      blockMeshDict;
    }}
    // * * 

    convertToMeters {convert_to_meters};

    vertices
    (
    """

    # Add each vertex in order
    for v in vertices:
        bm_content += f"    ( {v[0]:>7.4f}  {v[1]:>7.4f}  {v[2]:>7.4f} )\n"

    bm_content += textwrap.dedent(f"""
    );

    blocks
    (
        hex (0 1 2 3 4 5 6 7) ({block_size[0]} {block_size[1]} {block_size[2]}) simpleGrading ({simple_grading[0]} {simple_grading[1]} {simple_grading[2]})
    );

    edges
    (
    );

    patches
    (
        // no patches defined here; snappyHexMesh will carve the block
    );

    mergePatchPairs
    (
    );
    """)

    with open(f"{prefix}/mesh/system/blockMeshDict", "w") as f:
        f.write(textwrap.dedent(bm_content).lstrip())

    # ─────────────────────────────────────────────────────────────────────────────
    # 3) Write createPatchDict
    cp_content = f"""
    /*--------------------------------*- C++ -*----------------------------------*\\
    | =========                 |                                                 |
    | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\\\    /   O peration     | Version:  2412                                  |
    |   \\\\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\\\/     M anipulation  |                                                 |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        location    "system";
        object      createPatchDict;
    }}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    pointSync false;

    patches
    (
        {{
            name inlet;
            patchInfo
            {{
                type patch;
            }}
            constructFrom patches;
            patches (vesselSurface);
            filterFaces
            {{
                type box;
                box (({inlet_box[0][0]} {inlet_box[0][1]} {inlet_box[0][2]})  ({inlet_box[1][0]} {inlet_box[1][1]} {inlet_box[1][2]}));
            }}
        }}
        
        {{
            name outlet;
            patchInfo
            {{
                type patch;
            }}
            constructFrom patches;
            patches (vesselSurface);
            filterFaces
            {{
                type box;
                box (({outlet_box[0][0]} {outlet_box[0][1]} {outlet_box[0][2]})  ({outlet_box[1][0]} {outlet_box[1][1]} {outlet_box[1][2]}));
            }}
        }}
    );
    """

    with open(f"{prefix}/mesh/system/createPatchDict", "w") as f:
        f.write(textwrap.dedent(cp_content).lstrip())

    # ─────────────────────────────────────────────────────────────────────────────
    # 4) Write fvSolution (empty in your example)
    fv_solution = """
    /*--------------------------------*- C++ -*----------------------------------*\\
    | =========                 |                                                 |
    | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\\\    /   O peration     | Version:  2412                                  |
    |   \\\\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\\\/     M anipulation  |                                                 |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {
        version     2.0;
        format      ascii;
        class       dictionary;
        object      fvSolution;
    }
    // * * 

    solvers
    {
        // empty: no equation solves right now
    }

    relaxations
    {
        // empty
    }
    """
    with open(f"{prefix}/mesh/system/fvSolution", "w") as f:
        f.write(textwrap.dedent(fv_solution).lstrip())

    # ─────────────────────────────────────────────────────────────────────────────
    # 5) Write controlDict
    cd_content = f"""
    /*--------------------------------*- C++ -*----------------------------------*\\
    | =========                 |                                                 |
    | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\\\    /   O peration     | Version:  2412                                  |
    |   \\\\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\\\/     M anipulation  |                                                 |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      controlDict;
    }}
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    application     {application};

    startFrom       startTime;
    startTime       {startTime};

    stopAt          endTime;
    endTime         {endTime};

    deltaT          {deltaT};

    writeControl    timeStep;
    writeInterval   {writeInterval};

    purgeWrite      0;
    writeFormat     ascii;
    writePrecision  6;
    writeCompression off;

    timeFormat      general;
    timePrecision   6;

    runTimeModifiable true;
    """
    with open(f"{prefix}/mesh/system/controlDict", "w") as f:
        f.write(textwrap.dedent(cd_content).lstrip())

    # ─────────────────────────────────────────────────────────────────────────────
    # 6) Write fvSchemes
    fv_schemes = f"""
    /*--------------------------------*- C++ -*----------------------------------*\\
    | =========                 |                                                 |
    | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\\\    /   O peration     | Version:  2412                                  |
    |   \\\\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\\\/     M anipulation  |                                                 |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      fvSchemes;
    }}
    // * * 

    ddtSchemes
    {{
        default         {ddtScheme};
    }}

    gradSchemes
    {{
        default         {gradScheme};
    }}

    divSchemes
    {{
        default         {divScheme};
    }}

    laplacianSchemes
    {{
        default         {laplacianScheme};
    }}

    interpolationSchemes
    {{
        default         {interpolationScheme};
    }}

    snGradSchemes
    {{
        default         {snGradScheme};
    }}

    fluxRequired
    {{
        default         no;
    }}
    """
    with open(f"{prefix}/mesh/system/fvSchemes", "w") as f:
        f.write(textwrap.dedent(fv_schemes).lstrip())

    # ─────────────────────────────────────────────────────────────────────────────
    # 7) Write meshQualityDict
    mq_content = f"""
    /*--------------------------------*- C++ -*----------------------------------*\\
    | =========                 |                                                 |
    | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\\\    /   O peration     | Version:  2412                                  |
    |   \\\\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\\\/     M anipulation  |                                                 |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      meshQualityDict;
    }}
    // * * 

    maxNonOrtho           {maxNonOrtho};
    maxBoundarySkewness   {maxBoundarySkewness};
    maxInternalSkewness   {maxInternalSkewness};
    maxConcave            {maxConcave};
    minFlatness           {minFlatness};
    minVol                {minVol};
    minTetQuality         {minTetQuality};
    minArea               {minArea};
    minTwist              {minTwist};
    minDeterminant        {minDeterminant};
    minFaceWeight         {minFaceWeight};
    minVolRatio           {minVolRatio};
    minTriangleTwist      {minTriangleTwist};
    nSmoothScale          {nSmoothScale};
    errorReduction        {errorReduction};

    // (end of meshQualityDict)
    """
    with open(f"{prefix}/mesh/system/meshQualityDict", "w") as f:
        f.write(textwrap.dedent(mq_content).lstrip())

    # ─────────────────────────────────────────────────────────────────────────────
    # 8) Write snappyHexMeshDict
    shm_content = f"""
    /*--------------------------------*- C++ -*----------------------------------*\\
    | =========                 |                                                 |
    | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\\\    /   O peration     | Version:  2412                                  |
    |   \\\\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\\\/     M anipulation  |                                                 |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      snappyHexMeshDict;
    }}
    // * * 

    castellatedMesh {str(castellatedMesh).lower()};
    snap            {str(snap).lower()};
    addLayers       {str(addLayers).lower()};   // no prism layers for now

    geometry
    {{
        {stl_filename}
        {{
            type triSurfaceMesh;
            name vesselSurface;
        }}
    }}

    castellatedMeshControls
    {{
        maxLocalCells      {maxLocalCells};
        maxGlobalCells     {maxGlobalCells};
        minRefinementCells 0;

        // Allow up to {nCellsBetweenLevels} levels of refinement (coarse={refinement_levels[0]}, fine={refinement_levels[1]})
        nCellsBetweenLevels {nCellsBetweenLevels};

        features
        (
            // We did NOT run surfaceFeatureExtract, so leave empty
        );

        refinementSurfaces
        {{
            vesselSurface
            {{
                level ({refinement_levels[0]} {refinement_levels[1]});
            }}
        }};

        refinementRegions
        {{}};

        resolveFeatureAngle {resolveFeatureAngle}; // degrees
        allowFreeStandingZoneFaces false;
        locationInMesh ({locationInMesh[0]} {locationInMesh[1]} {locationInMesh[2]});
    }}

    snapControls
    {{
        nSmoothPatch 3;
        tolerance    1.0;   // mm (because convertToMeters = 0.001)
        nSolveIter   20;
        nRelaxIter   5;
    }}

    addLayersControls
    {{
        // no entries because addLayers = false
    }}

    meshQualityControls
    {{
        #include "meshQualityDict"
    }}

    debug            0;
    mergeTolerance   {mergeTolerance};
    """
    with open(f"{prefix}/mesh/system/snappyHexMeshDict", "w") as f:
        f.write(textwrap.dedent(shm_content).lstrip())

    # ─────────────────────────────────────────────────────────────────────────────
    # 9) Write surfaceFeatureExtractDict
    sfe_content = f"""
    /*--------------------------------*- C++ -*----------------------------------*\\
    | =========                 |                                                 |
    | \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
    |  \\\\    /   O peration     | Version:  2412                                  |
    |   \\\\  /    A nd           | Website:  www.openfoam.com                      |
    |    \\\\/     M anipulation  |                                                 |
    \\*---------------------------------------------------------------------------*/
    FoamFile
    {{
        version     2.0;
        format      ascii;
        class       dictionary;
        object      surfaceFeatureExtractDict;
    }}
    // * * 

    extractionMethod   extractFromSurface;
    extractFromSurfaceCoeffs
    {{
        // Name of input surface STL in constant/triSurface
        inputSurface        {input_surface};

        // Output edges will live in constant/triSurface as .eMesh
        writeObj            {str(writeObj).lower()};
        writeSTL            {str(writeSTL).lower()};
        writeThree          {str(writeThree).lower()};
        writeASCII          {str(writeASCII).lower()};

        triContinuity       {triContinuity};
        mergeTolerance      {mergeTolerance};
    }}
    """
    with open(f"{prefix}/mesh/system/surfaceFeatureExtractDict", "w") as f:
        f.write(textwrap.dedent(sfe_content).lstrip())

    # ─────────────────────────────────────────────────────────────────────────────
    print("All dictionary files created/overwritten in ./system/")
    print("blockMesh             # to build the background block")
    print("surfaceFeatureExtract # to extract edges (if needed)")
    print("snappyHexMesh         # to carve the block around your geometry")
    
    mesh_dir = os.path.join(prefix, "mesh")
    
    os.chdir(mesh_dir)
    ## print current directory
    print(f"Current directory: {os.getcwd()}")
    print(f"Running these commands in {mesh_dir}:\n")
    
    os.system(f". {FOAM_BASHRC} && blockMesh && surfaceFeatureExtract && snappyHexMesh -overwrite && foamToVTK")
    # os.system("blockMesh")
    # os.system("surfaceFeatureExtract")
    # os.system("snappyHexMesh -overwrite")
    # os.system("foamToVTK")
    # run_with_basic_runner(["blockMesh"])
    # run_with_basic_runner(["surfaceFeatureExtract"])
    # run_with_basic_runner(["snappyHexMesh", "-overwrite"])
    # run_with_basic_runner(["foamToVTK"])
    ## go back to the original directory
    print(f"Returning to original directory: {os.getcwd()}")
    os.chdir("../../")
    print("OpenFOAM mesh setup complete.")  

