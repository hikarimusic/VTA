import cv2
import numpy as np
from typing import List, Union

def merge_images_grid(image_grid: List[List[str]], output_path: str, font_size: float = 3) -> np.ndarray:
    """
    Merge multiple images in a grid layout based on a 2D array of image paths.
    Handles 'none' or empty slots in the grid.
    
    Parameters:
    image_grid (List[List[str]]): 2D list of image paths, with 'none' for empty slots
    output_path (str): Path where the merged image will be saved
    font_size (float): Size of the font for labels (default: 3)
    
    Returns:
    np.ndarray: The merged image array
    """
    # First pass: read all images and get dimensions
    processed_grid = []
    max_height_per_row = []
    max_width = 0
    
    for row_idx, row in enumerate(image_grid):
        processed_row = []
        row_max_height = 0
        row_total_width = 0
        
        for img_path in row:
            if img_path and img_path.lower() != 'none':
                img = cv2.imread(img_path)
                if img is None:
                    raise ValueError(f"Could not read image: {img_path}")
                height, width = img.shape[:2]
                row_max_height = max(row_max_height, height)
                row_total_width = max(row_total_width, width)
                processed_row.append(img)
            else:
                processed_row.append(None)
        
        processed_grid.append(processed_row)
        max_height_per_row.append(row_max_height)
        max_width = max(max_width, row_total_width)
    
    # Calculate total dimensions
    total_height = sum(max_height_per_row)
    total_width = max_width * max(len(row) for row in image_grid)
    
    # Create the output image
    merged_image = np.zeros((total_height, total_width, 3), dtype=np.uint8)
    merged_image.fill(255)  # Fill with white background
    
    # Place images and add labels
    current_y = 0
    labels = iter('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
    padding = 20
    font = cv2.FONT_HERSHEY_TRIPLEX
    thickness = 3
    
    for row_idx, (row, row_height) in enumerate(zip(processed_grid, max_height_per_row)):
        current_x = 0
        cell_width = total_width // len(row)
        
        for col_idx, img in enumerate(row):
            if img is not None:
                # Center the image in its cell
                h, w = img.shape[:2]
                y_offset = current_y + (row_height - h) // 2
                x_offset = current_x + (cell_width - w) // 2
                
                # Place the image
                merged_image[y_offset:y_offset + h, x_offset:x_offset + w] = img
                
                # Add label
                label = next(labels)
                cv2.putText(merged_image, label,
                           (x_offset + padding, y_offset + padding + 20 * font_size),
                           font, font_size, (0, 0, 0), thickness)
            
            current_x += cell_width
        
        current_y += row_height
    
    # Save the merged image
    cv2.imwrite(output_path, merged_image)
    return merged_image

image_grid = [
    ['up_HALLMARK_MITOTIC_SPINDLE.png', 'up_HALLMARK_UV_RESPONSE_DN.png', 'up_HALLMARK_PROTEIN_SECRETION.png', 'up_HALLMARK_TGF_BETA_SIGNALING.png'],
    ['down_HALLMARK_OXIDATIVE_PHOSPHORYLATION.png', 'down_HALLMARK_DNA_REPAIR.png', 'down_HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY.png', None]
]
merged = merge_images_grid(image_grid, 'GSEA_cluster1.png')

image_grid = [
    ['up_HALLMARK_MYC_TARGETS_V1.png', 'up_HALLMARK_E2F_TARGETS.png', 'up_HALLMARK_G2M_CHECKPOINT.png', None],
    ['down_HALLMARK_XENOBIOTIC_METABOLISM.png', 'down_HALLMARK_BILE_ACID_METABOLISM.png', 'down_HALLMARK_FATTY_ACID_METABOLISM.png', 'down_HALLMARK_PEROXISOME.png'],
    ['down_HALLMARK_HEME_METABOLISM.png', 'down_HALLMARK_ANDROGEN_RESPONSE.png', 'down_HALLMARK_UV_RESPONSE_DN.png', None]
]
merged = merge_images_grid(image_grid, 'GSEA_cluster2.png')

image_grid = [
    ['up_HALLMARK_XENOBIOTIC_METABOLISM.png', 'up_HALLMARK_BILE_ACID_METABOLISM.png', 'up_HALLMARK_OXIDATIVE_PHOSPHORYLATION.png', 'up_HALLMARK_FATTY_ACID_METABOLISM.png'],
    ['up_HALLMARK_ADIPOGENESIS.png', 'up_HALLMARK_PEROXISOME.png', 'up_HALLMARK_COAGULATION.png', 'up_HALLMARK_CHOLESTEROL_HOMEOSTASIS.png'],
    ['down_HALLMARK_G2M_CHECKPOINT.png', 'down_HALLMARK_MITOTIC_SPINDLE.png', 'down_HALLMARK_INFLAMMATORY_RESPONSE.png', 'down_HALLMARK_APICAL_JUNCTION.png'],
    ['down_HALLMARK_KRAS_SIGNALING_UP.png', None, None, None]
]
merged = merge_images_grid(image_grid, 'GSEA_cluster3.png')

