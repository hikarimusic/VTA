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

# Usage
image_grid = [
    ['cluster_heatmap_median_RNAseq.png', 'cluster_heatmap_median_ssGSEA.png'],
    ['cluster_heatmap_mean_RNAseq.png', 'cluster_heatmap_mean_ssGSEA.png'],
    ['cluster_heatmap_none_RNAseq.png', 'cluster_heatmap_none_ssGSEA.png']
]
merged = merge_images_grid(image_grid, 'compare_normalization.png')

