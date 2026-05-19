import numpy as np


def generate_complex_normal(mean, std, size):
    """Generates complex numbers with normal distribution.

    Args:
        mean (Float): Mean noise value.
        std (Float): Standard deviation noise value. 
        size (Tuple): Shape of array produced. 

    Returns:
        Array: Complex noise array. 
    """
    real = np.random.normal(mean, std, size)
    imag = np.random.normal(mean, std, size)
    return real + 1j * imag


def calculate_snr(img, vec_noise):
    """Calculate the Signal-to-Noise Ratio (SNR) of an input image.

    Information:
    https://onlinelibrary.wiley.com/doi/10.1002/mrm.28114
    http://eknygos.lsmuni.lt/springer/339/49-61.pdf

    Args:
        img (Array): Image.
        vec_noise (Array): Noise vector of voxels selected from image. 

    Returns:
        img_snr: SNR image. 
    """
    img_snr = (np.sqrt(2 - (np.pi/2)))*(abs(img) -
                                        np.mean(abs(vec_noise)))/np.std(abs(vec_noise))

    return img_snr


def calculate_rmse(image1, image2):
    """Calculate the Root Mean Square Error (RMSE) between two input images.

    Args:
        image1 (Array): Image 1.
        image2 (Array): Image 2.

    Returns:
        Float: RMSE value. 
    """
    diff_squared = (image1 - image2) ** 2
    mean_squared_diff = np.mean(diff_squared)
    rmse = np.sqrt(mean_squared_diff)

    return rmse


def ir(image, k=3):
    """Returns an image slice rotated by 90 degrees k times. Helpful for plotting reconstructions with plt.imshow().

    Args:
        image (Array): CPU array containing an image.
        k (int): number of times image is rotated about axis.

    Returns:
        Array: Rotated array.
    """
    image_out = np.rot90(image, k)

    return image_out


def ps(x, dx=0, dy=0):
    """Returns a circshifted image.

    Args:
        x (Array): 2D CPU array containing an image.
        dx (int): Number of voxels shifted in x.
        dy (int): Number of voxels shifted in y.

    Returns:
        Array: Rotated array.
    """

    return np.roll(x, (dx, dy), axis=(0, 1))
