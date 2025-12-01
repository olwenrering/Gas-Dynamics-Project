#!/usr/bin/env python3
# coding: utf-8

# Code d'Evolution Stellaire, Adaptatif et Modulaire for the 2020 decade.
# Copyright (c) 1997-2023 The Cesam2k20 authors
# SPDX-License-Identifier : GPL-3.0-or-later

import zipfile
import os



class ZipManager:
    def __init__( self, zip_filename, mode ):
        self.mode         = mode
        self.zip_filename = zip_filename

    def __enter__(self):
        self.zipf = zipfile.ZipFile(self.zip_filename, self.mode)
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        self.zipf.close()

    def zip_files( self, file_list ):
        """
        Compresses a list of files into a zip file.
        """
        with zipfile.ZipFile(self.zip_filename, 'w') as zipf:
            for file in file_list:
                if os.path.isfile(file):
                    zipf.write(file, os.path.basename(file))
                else:
                    print(f"The file {file} does not exist and cannot be added to the zip.")

    def list_files( self ):
        """
        Lists the files contained in the zip file.
        """
        file_list = []
        with zipfile.ZipFile(self.zip_filename, 'r') as zipf:
            file_list = zipf.namelist()
        return file_list

    def read_utf8file( self, filename ):
        """
        Returns the content of a specific file in the zip file.
        """
        with zipfile.ZipFile(self.zip_filename, 'r') as zipf:
            if filename in zipf.namelist():
                with zipf.open( filename, mode='r' ) as file:
                    return file.read().decode('utf-8')  # Assumes the content is UTF-8 text
            else:
                raise FileNotFoundError(f"The file {filename} does not exist in the zip.")

    def read_FortranBinaryfile( self, filename ):
        """
        Returns the content of a specific file in the zip file.
        """
        with zipfile.ZipFile(self.zip_filename, 'r') as zipf:
            if filename in zipf.namelist():
                file = FortranBinaryFile( agsm_r )
                with zipf.open( filename, mode='r' ) as file:
                    return file.read().decode('utf-8')  # Assumes the content is UTF-8 text
            else:
                raise FileNotFoundError(f"The file {filename} does not exist in the zip.")