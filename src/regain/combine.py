#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 University of Utah

import pandas as pd
import os

def combine_matrices(matrix1_path, matrix2_path, output_filename='combined_matrix.csv', delete_duplicates=False):
    output_path = os.path.join(os.getcwd(), output_filename)

    if not os.path.isfile(matrix1_path):
        raise FileNotFoundError(f" \033[91mError: The file '{matrix1_path}' does not exist\033[0m")
    if not os.path.isfile(matrix2_path):
        raise FileNotFoundError(f" \033[91mError: The file '{matrix2_path}' does not exist\033[0m")

    try:
        matrix1 = pd.read_csv(matrix1_path, index_col=0)
        matrix2 = pd.read_csv(matrix2_path, index_col=0)
    except pd.errors.EmptyDataError:
        raise ValueError(f" \033[91mError: One of the input files is empty: '{matrix1_path}' or '{matrix2_path}'\033[0m")
    except pd.errors.EmptyDataError:
        raise ValueError(f" \033[91mError: One of the input files could not be parsed: '{matrix1_path}' or '{matrix2_path}'. Check the file format\033[0m")

    matrix1.index.name = 'file'
    matrix2.index.name = 'file'

    combined_matrix = pd.concat([matrix1, matrix2], axis=1).fillna(0).astype(int)

    duplicate_columns = combined_matrix.columns[combined_matrix.columns.duplicated()].tolist()
    if duplicate_columns and not delete_duplicates:
        print(f"\n \033[91mWarning: Duplicate column headers found in combined matrix: {duplicate_columns}\033[0m\n")
        print(f" \033[93mDuplicate values need to be removed prior to Bayesian Network Structure Learning Analysis\n Why don't we do this for you? Well, we don't want to assume the duplicates were a mistake!\n If you do want ReGAIN to delete all duplicate values unconditionally, try adding the\n '--delete-duplicates' flag \033[0m\n")

    if delete_duplicates:
        combined_matrix = combined_matrix.loc[:, ~combined_matrix.columns.duplicated()]
        print("\n \033[95mDuplicate columns have been removed from combined matrix\033[0m\n")

    try:
        combined_matrix.to_csv(output_path)
    except PermissionError:
        raise PermissionError(f" \033[91mError: The file '{output_path}' could not be saved. Check your permissions\033[0m")

    print(f" \033[92mCombined matrix saved to {output_path}\033[0m \n")

def combine_metadata(metadata1_path, metadata2_path, output_filename='combined_metadata.csv', delete_duplicates=False):

    output_path = os.path.join(os.getcwd(), output_filename)

    if not os.path.isfile(metadata1_path):
        raise FileNotFoundError(f" \033[91mError: The file '{metadata1_path}' does not exist\033[0m")
    if not os.path.isfile(metadata2_path):
        raise FileNotFoundError(f" \033[91mError: The file '{metadata2_path}' does not exist\033[0m")

    try:
        metadata1 = pd.read_csv(metadata1_path)
        metadata2 = pd.read_csv(metadata2_path)
    except pd.errors.EmptyDataError:
        raise ValueError(f" \033[91m Error: One of the input files is empty: '{metadata1_path}' or '{metadata2_path}'\033[0m")
    except pd.errors.ParserError:
        raise ValueError(f" \033[91m Error: One of the input files could not be parsed: '{metadata1_path}' or '{metadata2_path}'. Check the file format\033[0m")

    metadata1[metadata1.columns[1]] = metadata1[metadata1.columns[1]].fillna('unspecified')
    metadata2[metadata2.columns[1]] = metadata2[metadata2.columns[1]].fillna('unspecified')

    for df in [metadata1, metadata2]:
        for col in df.columns[:2]:
            df[col] = df[col].astype(str)

    combined_metadata = pd.merge(metadata1, metadata2, left_on=list(metadata1.columns[:2]), right_on=list(metadata2.columns[:2]), how='outer')

    duplicate_values = combined_metadata[combined_metadata.duplicated(subset=combined_metadata.columns[0], keep=False)]
    if not duplicate_values.empty and not delete_duplicates:
        print(f" \033[91mWarning: Duplicate values found in column 1 of combined metadata:\n\n{duplicate_values}\033[0m \n")
        print(f" \033[93mDuplicate values need to be removed prior to Bayesian Network Structure Learning Analysis\n Why don't we do this for you? Well, we don't want to assume the duplicates are a mistake!\n If you do want ReGAIN to delete all duplicate values unconditionally, try adding the\n '--delete-duplicates' flag \033[0m\n")

    if delete_duplicates:
        combined_metadata.sort_values(by=[metadata1.columns[1]], inplace=True)
        combined_metadata.drop_duplicates(subset=combined_metadata.columns[0], keep='first', inplace=True)
        print(" \033[95mDuplicate rows have been removed from combined metadata\033[0m \n")

    if combined_metadata.empty:
        raise ValueError(" \033[91mError: Combined metadata resulted in an empty dataframe after merging. Check the input data\033[0m")

    try:
        combined_metadata.to_csv(output_path, index=False)
    except PermissionError:
        raise PermissionError(f" \033[91mError: The file '{output_path}' could not be saved. Check your permissions\033[0m")

    print(f" \033[92mCombined metadata saved to {output_path}\033[0m\n")

def run(args):
    try:
        combine_matrices(args.matrix1, args.matrix2, delete_duplicates=args.delete_duplicates)
        combine_metadata(args.metadata1, args.metadata2, delete_duplicates=args.delete_duplicates)
    except Exception as e:
        print(f" \033[91mAn unexpected error occurred: {str(e)}\033[0m")