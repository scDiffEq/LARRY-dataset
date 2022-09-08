# import packages: ------------------------------------------------------------
import licorice_font
import os
import scipy

# main module class: ----------------------------------------------------------
class _mtx_to_npz:
    def __init__(self, mtx_filepath):
        """
        1. read mtx
        2. convert to npz
        3. read npz

        - if you have `.npz`, do (3)
        - if you do not have `.npz`, do (1) and (2), followed by (3)
        """
        self.mtx = False
        self.npz_filepath = False
        self.mtx_filepath = mtx_filepath
        note = licorice_font.font_format("NOTE", ["BLUE"])
        self._msg = " - {} | reading and converting .mtx to .npz".format(note)

    def npz_filename(self):
        self.npz_filepath = self.mtx_filepath.replace("mtx.gz", "npz")

    def check_npz(self):
        if not self.npz_filepath:
            self.npz_filename()
        npz_exists = os.path.exists(self.npz_filepath)
        return npz_exists

    def read_mtx(self):
        if not self.check_npz():
            print(self._msg)
            if not self.mtx:
                self.mtx = scipy.io.mmread(self.mtx_filepath)
            else:
                self.save_npz()

    def save_npz(self):
        scipy.sparse.save_npz(self.npz_filepath, self.mtx)

    def read_npz(self):
        if not self.check_npz():
            self.read_mtx()
            self.save_npz()
        return scipy.sparse.load_npz(self.npz_filepath)


# main module-controlling function : ------------------------------------------
def _read_mtx_as_npz(mtx_path):
    mtx_npz = _mtx_to_npz(mtx_path)
    return mtx_npz.read_npz()