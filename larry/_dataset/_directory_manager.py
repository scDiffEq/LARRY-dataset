
# -- import packages: ----------------------------------------------------------
from pathlib import Path
from licorice_font import font_format

# -- Main class: ---------------------------------------------------------------
class DirectoryManager:
    def __init__(self, path: str, silent: bool = False):

        self._path = path
        self._dirs_created = []
        self._silent = silent

    @property
    def split_path(self):
        return [subdir for subdir in self._path.split("/") if subdir]    

    def notify(self, path):
            msg = font_format(f"mkdir: {path}", ['BLUE'])
            print(f" - [NOTE] | {msg}")

    def __mkdir__(self, path):
        path.mkdir()
        self._dirs_created.append(path)
        if not self._silent:
            self.notify(path)


    def __call__(self):

        target_dir = []

        for subdir in self.split_path:
            target_dir.append(subdir)
            subdir_ = Path("/".join(target_dir)).absolute()
            if not subdir_.exists():
                self.__mkdir__(subdir_)


# -- API-facing function: ------------------------------------------------------
def mkdir(path: str, silent: bool = False)->None:
    dir_manager = DirectoryManager(path = path)
    dir_manager()
