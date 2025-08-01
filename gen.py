import argparse
import numpy as np
import nibabel as nib
import meshio
from scipy import ndimage
from skimage import measure
from vmtk import vmtkscripts
from generate_foam_setup import mesh_generation
from totalsegmentator.python_api import totalsegmentator
import plotly.graph_objects as go
import os
from stl import mesh
from vtk import vtkXMLPolyDataReader
from skimage.morphology import skeletonize,medial_axis
from PyFoam.RunDictionary.SolutionDirectory import SolutionDirectory
import matplotlib.pyplot as plt
import plotly.graph_objects as go

def load_nifti(path):       
    nii = nib.load(path)
    data = nii.get_fdata().astype(np.int32)
    spacing = nii.header.get_zooms()
    return data, spacing

def get_terminal_points(label_nii_path):
    mask = nib.load(label_nii_path).get_fdata()
    mask= (mask==1)
    skeleton = skeletonize(mask)
    binary_structure = ndimage.generate_binary_structure(3, 3)
    neighbor_count = ndimage.convolve(skeleton.astype(int), binary_structure.astype(int), mode='constant')
    terminal_points = (skeleton > 0) & (neighbor_count == 2)
    
    # skel,distance = medial_axis(image=mask,return_distance=True)
    
    return terminal_points,skeleton

def segment_aorta(img_path, output_path="aorta_segmentation"):
    # Run TotalSegmentator for aorta
    totalsegmentator(img_path, output_path, roi_subset=["aorta" ])
    aorta_vol, _ = load_nifti(f"{output_path}/aorta.nii.gz")
    return aorta_vol == 1


def intersect_masks(mask1, mask2):
    return np.logical_and(mask1, mask2)


def output_stl(binary_mask, stl_path,voxel_spacing=(1.0, 1.0, 1.0)):
    

    verts, faces, normals, values = measure.marching_cubes(
        volume=binary_mask,
        level=0.5,
        spacing=voxel_spacing
    )

    # Write that mesh out as an STL for OpenFOAM
    # Convert faces to int64 and wrap as a meshio.Mesh
    mesh = meshio.Mesh( 
        points=verts,
        cells=[("triangle", faces.astype(np.int64))]
    )
    stl_filename = "coronary_surface.stl"
    meshio.write(stl_filename, mesh, file_format="stl")
    print(f"STL written for OpenFOAM: {stl_filename}")


def select_largest_components(binary_mask, num_keep=2):
    labeled, count = ndimage.label(binary_mask)
    sizes = ndimage.sum(binary_mask, labeled, range(count + 1))
    keep = np.argsort(sizes)[-num_keep:]
    refined = np.isin(labeled, keep)
    centroids = ndimage.center_of_mass(binary_mask, labeled, keep)
    return refined, centroids


def visualize_segmentation(binary_mask, aorta_mask, centroids, output_html="coronary_aorta_segmentation.html"):
    fig = go.Figure()
    # coronary mask volume
    fig.add_trace(go.Volume(
        x=np.arange(binary_mask.shape[0]),
        y=np.arange(binary_mask.shape[1]),
        z=np.arange(binary_mask.shape[2]),
        value=binary_mask.flatten(),
        isomin=0, isomax=1,
        opacity=0.1, surface_count=1
    ))
    # aorta mask volume
    fig.add_trace(go.Volume(
        x=np.arange(aorta_mask.shape[0]),
        y=np.arange(aorta_mask.shape[1]),
        z=np.arange(aorta_mask.shape[2]),
        value=aorta_mask.flatten(),
        isomin=0, isomax=1,
        opacity=0.1, surface_count=1,
        colorscale='Reds'
    ))
    # centroids
    for c in centroids:
        fig.add_trace(go.Scatter3d(
            x=[c[0]], y=[c[1]], z=[c[2]],
            mode='markers', marker=dict(size=5, color='blue')
        ))
    fig.update_layout(
        title="Coronary & Aorta Segmentation",
        scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z')
    )
    fig.write_html(output_html)
    print(f"Visualization saved to {output_html}")


def extract_surface_mesh(binary_mask, spacing):
    verts, faces, _, _ = measure.marching_cubes(
        volume=binary_mask, level=0.5, spacing=spacing
    )
    mesh = meshio.Mesh(points=verts, cells=[("triangle", faces.astype(np.int64))])
    return mesh, verts


def write_stl(mesh, stl_path):
    meshio.write(stl_path, mesh, file_format="stl")
    print(f"STL written to {stl_path}")


def bounding_box(stl_path:str):
    stl_mesh = mesh.Mesh.from_file(stl_path)
    min_coords = np.min(stl_mesh.points.reshape(-1, 3), axis=0)
    max_coords = np.max(stl_mesh.points.reshape(-1, 3), axis=0)
    return min_coords, max_coords

def extract_min_max(verts):
    return np.min(verts, axis=0), np.max(verts, axis=0)


def generate_openfoam_setup(stl_path, min_coords, max_coords,prefix=""):
    vertices = np.array([
        [min_coords[0], min_coords[1], min_coords[2]],
        [max_coords[0], min_coords[1], min_coords[2]],
        [max_coords[0], max_coords[1], min_coords[2]],
        [min_coords[0], max_coords[1], min_coords[2]],
        [min_coords[0], min_coords[1], max_coords[2]],
        [max_coords[0], min_coords[1], max_coords[2]],
        [max_coords[0], max_coords[1], max_coords[2]],
        [min_coords[0], max_coords[1], max_coords[2]]
    ])
    mesh_generation(stl_filename=stl_path, vertices=vertices,prefix=prefix)




# METHOD 4: Interactive HTML with terminal points highlighted
def save_skeleton_html_with_terminals(skeleton,terminals, output_file="skeleton_with_terminals.html"):    
    skel_coords = np.where(skeleton > 0)
    term_coords = np.where(terminals > 0)    
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=skel_coords[0],
        y=skel_coords[1],
        z=skel_coords[2],
        mode='markers',
        marker=dict(size=2, color='red', opacity=0.6),
        name='Skeleton'
    ))
    
    # Add terminal points
    fig.add_trace(go.Scatter3d(
        x=term_coords[0],
        y=term_coords[1],
        z=term_coords[2],
        mode='markers',
        marker=dict(size=5, color='blue', opacity=0.9),
        name='Terminal Points'
    ))
    fig.update_layout(title='Skeleton with Terminal Points')
    fig.write_html(output_file)
    print(f"Skeleton with terminals saved to {output_file}")
    print(f"Found {len(term_coords[0])} terminal points")

def read_vtp(filename):
    reader = vtkXMLPolyDataReader()
    reader.SetFileName(filename)
    reader.Update()
    return reader.GetOutput()


        

def main(args):
    
    ## making foam directories
    os.makedirs(f"{args.prefix}/mesh/system", exist_ok=True)
    os.makedirs(f"{args.prefix}/mesh/constant", exist_ok=True)
    os.makedirs(f"{args.prefix}/mesh/constant/polySurface", exist_ok=True)
    os.makedirs(f"{args.prefix}/mesh/constant/triSurface", exist_ok=True)
    os.makedirs(f"{args.prefix}/mesh/0", exist_ok=True)
    
    # Load volumes
    label_vol, spacing = load_nifti(args.label_nii_path)
    coronary_mask = (label_vol == 1)
    output_stl(coronary_mask, args.stl_output_path, voxel_spacing=spacing)

    terminal_points,skeleton = get_terminal_points(args.label_nii_path)
    # save_skeleton_html_with_terminals(skeleton, terminal_points, output_file=f"{args.prefix}/skeleton_with_terminals.html")
    
    print(f"Found {np.sum(terminal_points[0])} terminal points in the coronary mask.")
    
    
    
    # Aorta segmentation & intersection
    aorta_mask = segment_aorta(args.img_nii_path,output_path=f"{args.prefix}/aorta_segmentation")
    combined_mask = intersect_masks(coronary_mask, aorta_mask)
    ## making a plotly mesh with the combined mask in red, and the coronary mask in blue
    visualize_segmentation(combined_mask, aorta_mask, [], output_html=f"{args.prefix}/coronary_aorta_segmentation.html")
    os.system(f"cp {args.stl_output_path} {args.prefix}/mesh/constant/triSurface/coronary_surface.stl") 
    min_c, max_c = bounding_box(args.stl_output_path)
    print(f"Bounding box min: {min_c}, max: {max_c}")
    generate_openfoam_setup(args.stl_output_path, min_c, max_c,prefix=args.prefix)

    ## converting the terminal points into mesh coordinates.
    mesh = read_vtp(f"{args.prefix}/mesh/VTK/mesh_0/boundary/vesselSurface.vtp")
    print(mesh.GetOrigin())
    print(type(mesh))
    terminal_points_on_mesh = []
    for i, j, k in np.argwhere(terminal_points[0]):
        x = i * spacing[0] + mesh.GetOrigin()[0]
        y = j * spacing[1] + mesh.GetOrigin()[1]
        z = k * spacing[2] + mesh.GetOrigin()[2]
        pid = mesh.FindPoint((x, y, z))
        terminal_points_on_mesh.append(mesh.GetPoint(pid))
    terminal_points_on_mesh = np.array(terminal_points_on_mesh)
    print(f"Terminal points on mesh: {terminal_points_on_mesh.shape[0]} points.")
    
    ## creating a clone of the mesh directory, but now for solving 
    foam_dir = SolutionDirectory(os.path.join(args.prefix,"mesh"), archive=None)
    foam_dir.clone(os.path.join(args.prefix,"foam"), archive=None)
    # Write terminal points to a file
    terminal_points_file = os.path.join(args.prefix, "foam", "constant", "terminal_points.txt")
    np.savetxt(terminal_points_file, terminal_points_on_mesh, fmt='%.6f', header='x y z')
    print(f"Terminal points saved to {terminal_points_file}")
 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Coronary mesh & simulation setup pipeline"
    )
    parser.add_argument("--label_nii_path", required=True)
    parser.add_argument("--img_nii_path", required=True)
    parser.add_argument(
        "--stl_output_path", default="coronary_surface.stl"
    )
    parser.add_argument(
        "--prefix",type=str, 
        required=True,
    )
    args = parser.parse_args()
    main(args)
