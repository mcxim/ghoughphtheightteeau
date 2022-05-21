import os.path
from binary_file_search.BinaryFileSearch import BinaryFileSearch
import re
from cytoolz import memoize
from string import ascii_lowercase
from collections import defaultdict
import sys

PREFER_BIGGER_STEPS = 1.0

FORCE_END_MATCH = True
FORCE_START_MATCH = True


def lowercase_only(string):
    return "".join(char for char in string if char in ascii_lowercase)


def parse_table():
    table = dict()
    with open(
        os.path.join(
            os.path.dirname(os.path.realpath(__file__)),
            "cmudict.txt.m-mAlign.2-2.delX.1-best.conYX.align",
        ),
        "r",
    ) as fd:
        for line_ in fd:
            line = line_.strip()
            if line:
                [unparsed_graphemes, unparsed_phonemes] = line.strip().split("\t")
                table[unparsed_phonemes] = unparsed_graphemes
    return table


def ghoughphtheightteeau(word, max_len=float("inf")):
    result, as_ins = perform_ghoughphtheightteeau(
        tuple(get_arpabet(word)), max_len, start=True
    )
    print(result, end="")
    print("".join("\n\t{}".format(as_in) for as_in in as_ins))


def find_bars(string):
    return [m.start() for m in re.finditer(r"\|", string)]


def get_example(phonemes_option, graphemes_option, match):
    grapheme_idxs = [0] + [bar_idx + 1 for bar_idx in find_bars(graphemes_option)]
    phoneme_idxs = [0] + [bar_idx + 1 for bar_idx in find_bars(phonemes_option)]
    subword_start_idx = next(idx for idx in phoneme_idxs if idx >= match.start())
    subword_end_idx = next(idx for idx in reversed(phoneme_idxs) if idx <= match.end())
    return graphemes_option[
        grapheme_idxs[phoneme_idxs.index(subword_start_idx)] : grapheme_idxs[
            phoneme_idxs.index(subword_end_idx)
        ]
    ]


@memoize
def perform_ghoughphtheightteeau(phonemes, max_len, table=parse_table(), start=False):
    if not len(phonemes):
        return ("", tuple())
    options = []
    for step_len in range(1, min(max_len, len(phonemes)) + 1):
        pattern = r"(^|\|)" + r".(_\|)*".join(phonemes[:step_len]) + r"\|"
        if FORCE_END_MATCH:
            if step_len == len(phonemes):
                pattern += r"(?=(_\|)*$)"
            else:
                pattern += r"(?=[A-Z_])"
        if FORCE_START_MATCH and start:
            pattern = r"^" + pattern
        best_options, best_len = defaultdict(int), 0
        as_in = dict()
        for phonemes_option in table.keys():
            for match in re.finditer(pattern, phonemes_option):
                graphemes_option = table[phonemes_option]
                example_ = get_example(phonemes_option, graphemes_option, match)
                example = lowercase_only(example_)
                if len(example) >= best_len:
                    if len(example) > best_len:
                        best_options.clear()
                    best_options[example] += 1
                    best_len = len(example)
                    if example not in as_in:
                        as_in[example] = lowercase_only(graphemes_option)
        if not best_options:
            continue
        best_option = sorted(
            best_options.items(), key=lambda item: item[1], reverse=True
        )[0][0]
        best_rest, rest_as_ins = perform_ghoughphtheightteeau(
            phonemes[step_len:], max_len
        )
        options.append(
            (
                len(best_option + best_rest) + PREFER_BIGGER_STEPS * step_len,
                best_option + best_rest,
                "'{}' as in '{}'".format(best_option, as_in[best_option]),
                rest_as_ins,
            )
        )
    (score, best_overall_option, best_overall_as_ins, best_rest_as_ins) = sorted(
        options, key=lambda option: option[0], reverse=True
    )[0]
    return (best_overall_option, (best_overall_as_ins,) + best_rest_as_ins)


def get_arpabet(word):
    with BinaryFileSearch("./cmudict.txt", sep="\t", string_mode=True) as bfs:
        return bfs.search(" ".join(list(word)))[0][1].split(" ")


ghoughphtheightteeau(sys.argv[1])
