#!/bin/python3
import argparse
import os
import string
from enum import Enum


def doInline(filepath, seen=set()):
    if filepath in seen:
        return "", seen
    seen.add(filepath)
    
    outcode = f"// BEGIN PACKED FROM {filepath}\n"

    for line in open(filepath, "r").readlines():
        if line.startswith("#include \""):
            libraryname = line.partition("\"")[2].partition("\"")[0]
            filebase, _ = os.path.split(filepath)
            librarypath = os.path.join(filebase, libraryname)
            currcode, seen = doInline(librarypath, seen)
            if "\ufeff" in currcode:
                print(currcode)
                raise Exception
            outcode += currcode

        elif not line.startswith("#pragma once"):
            outcode += line
            if "\ufeff" in line:
                print(line)
                raise Exception
        

    outcode += "\n"
    outcode += f"// END PACKED FROM {filepath}\n"
    outcode += "\n"
    return outcode, seen


def lexer(code):
    code = code.replace("\t", "")
    for _ in range(30):
        code = code.replace("  ", " ")

    match = {"//": "\n", "/*": "*/", "\"": "\"", "'": "'", "#": "\n"}
    mode = "none"

    lexemes = [""]
    while code:
        if mode == "none":
            starts = [k for k in match if code.startswith(k)]
            if starts:
                assert len(starts) == 1
                mode = starts[0]
                lexemes.append(mode)
                code = code[len(mode):]
            else:
                lexemes[-1] += code[0]
                code = code[1:]
        else:
            if code.startswith(match[mode]):
                lexemes[-1] += match[mode]
                lexemes.append("")
                code = code[len(match[mode]):]
                mode = "none"
            else:
                lexemes[-1] += code[0]
                code = code[1:]
    lexemes = [x.replace("\n", "") for x in lexemes]
    lexemes = [x for x in lexemes if x not in ["", " "]]
    return lexemes


def remComments(code):
    lexemes = [x for x in lexer(code) if not x.startswith(
        "//") and not x.startswith("/*")]
    out = ["\n"+x+"\n" if x.startswith("#") else x for x in lexemes]
    return "".join(out).replace("\n\n", "\n").strip("\n")


def remWhitespace(code):
    varchars = string.ascii_letters + string.digits + "_" + "'\""
    out = ""
    for lexeme in lexer(code):
        line = ""
        for prev, curr, after in zip(lexeme, lexeme[1:], lexeme[2:]):
            if curr == " " and (prev not in varchars or after not in varchars):
                pass
            else:
                line += curr
        if lexeme[0] == "#":
            out += "\n" + lexeme[0] + line + lexeme[-1] + "\n"
        else:
            out += lexeme[0] + line + lexeme[-1]
    return out.replace("\n\n", "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="inline all custom library imports"
    )
    parser.add_argument(
        "filename",
        metavar="filename",
        type=str,
        help="the file path of the c++ code to be packed"
    )
    parser.add_argument(
        "-o", dest="output",
        type=str,
        default="packed.cpp",
        help="destination file for pack"
    )
    parser.add_argument(
        "--no_comments",
        action="store_true"
    )
    parser.add_argument(
        "--no_whitespace",
        action="store_true"
    )
    args = parser.parse_args()

    filepath = os.path.join(os.getcwd(), args.filename)
    code, libraries = doInline(filepath)
    #print(libraries)
    if args.no_comments:
        code = remComments(code)
    if args.no_whitespace:
        code = remWhitespace(code)
        
    assert args.output != args.filename
    open(args.output, "w").write(code)
