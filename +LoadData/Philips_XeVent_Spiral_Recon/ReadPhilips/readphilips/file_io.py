# -*- coding: utf-8 -*-
import tkinter as tk
from tkinter import filedialog

class io:
    """This module contains an abstract class for selecting files and folders.
    """

    def __init__(self, initial_dir="C:\\Users", ext_type="*.*", f_type = "All files"):
        self.initial_dir = initial_dir
        self.ext_type = ext_type
        self.f_type = f_type
        

    def selectFile(self):
        """Creates a file-explorer box to select a file.
        
        Usage:
            from readphilips.file_io import io
            io(f_type = "sin files", ext_type = "*.sin").selectFile()
        
        Args:
            initial_dir (string): Initial directory.
            ext_type (string): Extension (e.g. "*.txt".
            f_type (string): File type.
        
        Returns:
            filepath (string): Selected file name and path.
        """
        root = tk.Tk()
        root.withdraw()
        
        self.filepath = filedialog.askopenfilename(initialdir = self.initial_dir,
                                              title = "Select file",
                                              filetypes= [(self.f_type, self.ext_type)])
        return self.filepath
    
    
    def selectDirectory(self):
        """Creates a file-explorer box to select a folder.
        
        Usage:
            from readphilips.file_io import io
            io(initial_dir="C:\\Users").selectDirectory()
        
        Args:
            initial_dir (string): Initial directory.
                    
        Returns:
            directory (string): Selected file name and path.
        """
        root = tk.Tk()
        root.withdraw()
        self.directory = filedialog.askdirectory(initialdir = self.initial_dir,
                                                 title = "Select folder")
        
        return self.directory
    
