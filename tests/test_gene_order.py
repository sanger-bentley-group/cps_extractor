import pytest
import pandas as pd

from lib.check_gene_order import CheckGeneOrder


@pytest.fixture
def gene_checker():
    gene_check = CheckGeneOrder("tests/test_data/18B_test.gff3", "18B.gff")
    return gene_check


@pytest.fixture
def comparison_df():
    comparison_df = pd.read_csv("tests/test_data/df.csv", keep_default_na=False)
    return comparison_df


@pytest.fixture
def data_df():
    data_df = pd.read_csv("tests/test_data/df.csv", keep_default_na=False)
    data_df.drop("equal_gene", axis=1, inplace=True)
    return data_df


def test_sort_for_unifref100(gene_checker):
    uniref_sorted = gene_checker.sort_for_uniref("UniRef:UniRef100_A0A9Q9R4P9")
    assert uniref_sorted == (0, "UniRef:UniRef100_A0A9Q9R4P9")


def test_sort_for_unifref100_not_uniref(gene_checker):
    uniref_sorted = gene_checker.sort_for_uniref("blah")
    assert uniref_sorted == (1, "blah")


def test_get_features_not_ref(gene_checker):
    features = gene_checker.get_features("tests/test_data/truncated.gff3", False)
    assert features == [
        {
            "location": "[0:1446](+)",
            "name": "cps2a",
            "product": "Anionic cell wall polymer biosynthesis enzyme TagV/TagU, LytR-Cps2A-Psr (LCP) family (peptidoglycan teichoic acid transferase)",
            "uniref_id": "UniRef100_Q9AH98",
        }
    ]


def test_get_features_ref(gene_checker):
    features = gene_checker.get_features("tests/test_data/truncated.gff3", True)
    assert features == [
        {
            "location_ref": "[0:1446](+)",
            "name_ref": "cps2a",
            "product_ref": "Anionic cell wall polymer biosynthesis enzyme TagV/TagU, LytR-Cps2A-Psr (LCP) family (peptidoglycan teichoic acid transferase)",
            "uniref_id_ref": "UniRef100_Q9AH98",
        }
    ]


def test_get_features_no_gene_info(gene_checker):
    features = gene_checker.get_features(
        "tests/test_data/truncated_no_gene.gff3", False
    )
    assert features == [
        {
            "location": "[0:1446](+)",
            "name": "unidentified",
            "product": "Anionic cell wall polymer biosynthesis enzyme TagV/TagU, LytR-Cps2A-Psr (LCP) family (peptidoglycan teichoic acid transferase)",
            "uniref_id": "NA",
        }
    ]


def test_combine_ref_sample_data(gene_checker):
    combined_data = gene_checker.combine_ref_sample_data(
        [
            {
                "location": "[0:1446](+)",
                "name": "unidentified",
                "product": "blah",
                "uniref_id": "NA",
            }
        ],
        [
            {
                "location_ref": "[0:1446](+)",
                "name_ref": "cps2a",
                "product_ref": "blah",
                "uniref_id_ref": "UniRef100_Q9AH98",
            }
        ],
    )

    assert combined_data == [
        {
            "location": "[0:1446](+)",
            "name": "unidentified",
            "product": "blah",
            "uniref_id": "NA",
            "location_ref": "[0:1446](+)",
            "name_ref": "cps2a",
            "product_ref": "blah",
            "uniref_id_ref": "UniRef100_Q9AH98",
        }
    ]


def test_combine_ref_sample_data_extra_ref_gene(gene_checker):
    combined_data = gene_checker.combine_ref_sample_data(
        [
            {
                "location": "[0:1446](+)",
                "name": "unidentified",
                "product": "blah",
                "uniref_id": "NA",
            }
        ],
        [
            {
                "location_ref": "[0:1446](+)",
                "name_ref": "cps2a",
                "product_ref": "blah",
                "uniref_id_ref": "UniRef100_Q9AH98",
            },
            {
                "location_ref": "[0:1446](+)",
                "name_ref": "cps2a",
                "product_ref": "blah",
                "uniref_id_ref": "UniRef100_Q9AH98",
            },
        ],
    )

    assert combined_data == [
        {
            "location_ref": "[0:1446](+)",
            "name_ref": "cps2a",
            "product_ref": "blah",
            "uniref_id_ref": "UniRef100_Q9AH98",
            "location": "[0:1446](+)",
            "name": "unidentified",
            "product": "blah",
            "uniref_id": "NA",
        },
        {
            "location_ref": "[0:1446](+)",
            "name_ref": "cps2a",
            "product_ref": "blah",
            "uniref_id_ref": "UniRef100_Q9AH98",
        },
    ]


def test_check_gene_in_ref_and_sample_equal(gene_checker, comparison_df):
    row_in_df = [
        "glf",
        "UniRef100_Q4K161",
        "UDP-galactopyranose mutase",
        "[19746:20250](+)",
        "glf",
        "UniRef100_Q4K161",
        "UDP-galactopyranose mutase",
        "[19748:20252](+)",
        "Y",
    ]
    row_series = pd.Series(
        row_in_df,
        index=[
            "name",
            "uniref_id",
            "product",
            "location",
            "name_ref",
            "uniref_id_ref",
            "product_ref",
            "location_ref",
            "equal_gene",
        ],
    )
    assert gene_checker.check_gene_in_ref_and_sample(row_series, comparison_df) == "Y"


def test_check_gene_in_ref_na(gene_checker, comparison_df):
    row_in_df = [
        "glf",
        "NA",
        "UDP-galactopyranose mutase",
        "[19746:20250](+)",
        "glf",
        "UniRef100_Q4K161",
        "UDP-galactopyranose mutase",
        "[19748:20252](+)",
        "Y",
    ]
    row_series = pd.Series(
        row_in_df,
        index=[
            "name",
            "uniref_id",
            "product",
            "location",
            "name_ref",
            "uniref_id_ref",
            "product_ref",
            "location_ref",
            "equal_gene",
        ],
    )
    assert gene_checker.check_gene_in_ref_and_sample(row_series, comparison_df) == "N"


def test_check_gene_in_ref_and_sample_not(gene_checker, comparison_df):
    row_not_in_df = [
        "unidentified",
        "UniRef100_UPI0005E19987",
        "acyltransferase family protein",
        "[12213:12687](+)",
        "unidentified",
        "NA",
        "hypothetical protein",
        "[12213:12384](+)",
        "N",
    ]
    row_not_in_series = pd.Series(
        row_not_in_df,
        index=[
            "name",
            "uniref_id",
            "product",
            "location",
            "name_ref",
            "uniref_id_ref",
            "product_ref",
            "location_ref",
            "equal_gene",
        ],
    )
    assert (
        gene_checker.check_gene_in_ref_and_sample(row_not_in_series, comparison_df)
        == "N"
    )


def test_check_gene_in_ref_and_sample_dif_id(gene_checker, comparison_df):
    row_in_df = [
        "glf",
        "UniRef100_A0A9Q9R4P9",
        "UDP-galactopyranose mutase",
        "[19746:20250](+)",
        "glf",
        "UniRef100_Q4K161",
        "UDP-galactopyranose mutase",
        "[19748:20252](+)",
        "Y",
    ]
    row_in_series = pd.Series(
        row_in_df,
        index=[
            "name",
            "uniref_id",
            "product",
            "location",
            "name_ref",
            "uniref_id_ref",
            "product_ref",
            "location_ref",
            "equal_gene",
        ],
    )
    assert (
        gene_checker.check_gene_in_ref_and_sample(row_in_series, comparison_df) == "Y"
    )


def test_check_gene_in_ref_and_sample_not_name(gene_checker, comparison_df):
    row_not_in_df = [
        "blah",
        "blah",
        "UDP-galactopyranose mutase",
        "[19746:20250](+)",
        "glf",
        "UniRef100_Q4K161",
        "UDP-galactopyranose mutase",
        "[19748:20252](+)",
        "Y",
    ]
    row_not_in_series = pd.Series(
        row_not_in_df,
        index=[
            "name",
            "uniref_id",
            "product",
            "location",
            "name_ref",
            "uniref_id_ref",
            "product_ref",
            "location_ref",
            "equal_gene",
        ],
    )
    assert (
        gene_checker.check_gene_in_ref_and_sample(row_not_in_series, comparison_df)
        == "N"
    )


def test_gene_comparison(gene_checker, data_df):
    first_row_df = data_df.head(1)
    compared_df = gene_checker.gene_comparison(first_row_df)
    assert compared_df["equal_gene"].iloc[0] == "Y"
    assert compared_df["gene_in_ref_and_sample"].iloc[0] == "Y"


def test_gene_comparison_full(gene_checker, data_df):
    compared_df = gene_checker.gene_comparison(data_df)
    assert compared_df["equal_gene"].iloc[12] == "N"
    assert compared_df["gene_in_ref_and_sample"].iloc[12] == "N"
