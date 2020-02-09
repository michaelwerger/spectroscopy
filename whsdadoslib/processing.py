class ProcessingData(object):

    fits_path = None
    corr_path = None
    products_path = None

    @staticmethod
    def set_fits_path(fits_path):
        ProcessingData.fits_path = fits_path

    @staticmethod
    def set_corr_path(corr_path):
        ProcessingData.corr_path = corr_path

    @staticmethod
    def set_products_path(products_path):
        ProcessingData.products_path = products_path
    