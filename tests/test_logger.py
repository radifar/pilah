from pilah.logger import sort_records

sample_record = {
    ("A", 6, "GLU"): [4.18, {}, "Neutral"],
    ("A", 7, "GLU"): [4.56, {}, "Neutral"],
    ("A", 11, "GLU"): [4.7, {}, "Neutral"],
    ("A", 13, "LYS"): [11.08, {}, "Neutral"],
    ("A", 14, "GLU"): [4.15, {}, "Neutral"],
    ("A", 20, "ASP"): [2.89, {}, "Neutral"],
    ("A", 21, "LYS"): [10.77, {}, "Neutral"],
    ("A", 22, "ASP"): [4.32, {}, "Neutral"],
    ("A", 24, "ASP"): [4.53, {}, "Neutral"],
    ("A", 30, "LYS"): [10.52, {}, "Neutral"],
    ("A", 31, "GLU"): [4.46, {}, "Neutral"],
    ("A", 45, "GLU"): [3.8, {}, "Neutral"],
    ("A", 47, "GLU"): [3.48, {}, "Neutral"],
    ("A", 50, "ASP"): [3.71, {}, "Neutral"],
    ("A", 54, "GLU"): [3.7, {}, "Neutral"],
    ("A", 56, "ASP"): [3.67, {}, "Neutral"],
    ("A", 58, "ASP"): [4.16, {}, "Neutral"],
    ("A", 64, "ASP"): [3.63, {}, "Neutral"],
    ("A", 67, "GLU"): [5.4, {}, "Neutral"],
    ("A", 75, "LYS"): [11.17, {}, "Neutral"],
    ("A", 37, "ARG"): [14.0, {}, "Positive"],
    ("A", 74, "ARG"): [14.0, {}, "Positive"],
    ("A", 100, " NA"): [14.0, {}, "Positive"],
    ("A", 101, " MG"): [14.0, {}, "Positive"],
    ("A", 102, " ZN"): [14.0, {}, "Positive"],
    ("A", 103, " ZN"): [14.0, {}, "Positive"],
}

expected_sorted_sample_record = {
    ("A", 6, "GLU"): [4.18, {}, "Neutral"],
    ("A", 7, "GLU"): [4.56, {}, "Neutral"],
    ("A", 11, "GLU"): [4.7, {}, "Neutral"],
    ("A", 13, "LYS"): [11.08, {}, "Neutral"],
    ("A", 14, "GLU"): [4.15, {}, "Neutral"],
    ("A", 20, "ASP"): [2.89, {}, "Neutral"],
    ("A", 21, "LYS"): [10.77, {}, "Neutral"],
    ("A", 22, "ASP"): [4.32, {}, "Neutral"],
    ("A", 24, "ASP"): [4.53, {}, "Neutral"],
    ("A", 30, "LYS"): [10.52, {}, "Neutral"],
    ("A", 31, "GLU"): [4.46, {}, "Neutral"],
    ("A", 37, "ARG"): [14.0, {}, "Positive"],
    ("A", 45, "GLU"): [3.8, {}, "Neutral"],
    ("A", 47, "GLU"): [3.48, {}, "Neutral"],
    ("A", 50, "ASP"): [3.71, {}, "Neutral"],
    ("A", 54, "GLU"): [3.7, {}, "Neutral"],
    ("A", 56, "ASP"): [3.67, {}, "Neutral"],
    ("A", 58, "ASP"): [4.16, {}, "Neutral"],
    ("A", 64, "ASP"): [3.63, {}, "Neutral"],
    ("A", 67, "GLU"): [5.4, {}, "Neutral"],
    ("A", 74, "ARG"): [14.0, {}, "Positive"],
    ("A", 75, "LYS"): [11.17, {}, "Neutral"],
    ("A", 100, " NA"): [14.0, {}, "Positive"],
    ("A", 101, " MG"): [14.0, {}, "Positive"],
    ("A", 102, " ZN"): [14.0, {}, "Positive"],
    ("A", 103, " ZN"): [14.0, {}, "Positive"],
}


def test_sort_records():
    sorted_sample_record = sort_records(sample_record)

    assert sorted_sample_record == expected_sorted_sample_record
