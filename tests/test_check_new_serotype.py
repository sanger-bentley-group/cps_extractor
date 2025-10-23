import pytest
import pandas as pd

from lib.check_new_serotype import NewSerotype


@pytest.fixture
def new_serotype_no_disruptions():
    new_serotype = NewSerotype(
        "11F",
        "tests/test_data/known_disruptions.csv",
        "tests/test_data/disrupted_genes.csv",
    )
    return new_serotype


@pytest.fixture
def new_serotype_disruptions():
    new_serotype = NewSerotype(
        "4",
        "tests/test_data/known_disruptions.csv",
        "tests/test_data/disrupted_genes.csv",
    )
    return new_serotype


@pytest.fixture
def disrupted_gene_df(new_serotype_disruptions):
    disrupted_df = new_serotype_disruptions.read_csv(
        new_serotype_disruptions.disrupted_genes
    )
    return disrupted_df


@pytest.fixture
def known_disruptions_df(new_serotype_disruptions):
    known_disruption_df = new_serotype_disruptions.read_csv(
        new_serotype_disruptions.known_disruptions
    )
    return known_disruption_df


def test_read_csv(new_serotype_disruptions):
    data = {
        "sample": ["11511_7#17", "11511_7#16", "11511_7#15", "11511_7#15"],
        "gene": ["wcjE", "gct", "wcjE", "gct"],
        "gene_integrity": ["disrupted", "disrupted", "disrupted", "disrupted"],
    }
    actual = pd.DataFrame(data)
    df_read = new_serotype_disruptions.read_csv(
        new_serotype_disruptions.disrupted_genes
    )
    pd.testing.assert_frame_equal(actual, df_read)


def test_remove_known_disruptions_no_disruptions(
    new_serotype_no_disruptions, disrupted_gene_df, known_disruptions_df
):
    empty_df = new_serotype_no_disruptions.remove_known_disruptions(
        disrupted_gene_df, known_disruptions_df
    )
    assert empty_df.empty


def test_remove_known_disruptions(
    new_serotype_disruptions, disrupted_gene_df, known_disruptions_df
):
    data = {
        "sample": ["11511_7#17", "11511_7#16", "11511_7#15", "11511_7#15"],
        "gene": ["wcjE", "gct", "wcjE", "gct"],
        "gene_integrity": ["disrupted", "disrupted", "disrupted", "disrupted"],
    }
    actual = pd.DataFrame(data)
    filtered_df = new_serotype_disruptions.remove_known_disruptions(
        disrupted_gene_df, known_disruptions_df
    )
    pd.testing.assert_frame_equal(actual, filtered_df)


def test_find_unique_samples(new_serotype_disruptions, disrupted_gene_df):
    unique_samples = new_serotype_disruptions.find_unique_samples(disrupted_gene_df)
    assert unique_samples == ["11511_7#15", "11511_7#16", "11511_7#17"]


def test_find_unique_samples_non_unique(new_serotype_disruptions):
    data = {
        "sample": ["11511_7#17", "11511_7#16", "11511_7#15", "11511_7#15"],
        "gene": ["wcjE", "wcjE", "wcjE", "gct"],
        "gene_integrity": ["disrupted", "disrupted", "disrupted", "disrupted"],
    }
    non_unique_disruptions = pd.DataFrame(data)
    unique_samples = new_serotype_disruptions.find_unique_samples(
        non_unique_disruptions
    )
    assert unique_samples == ["11511_7#15", "11511_7#16"]


def test_find_unique_samples_no_disruptions(
    new_serotype_no_disruptions, disrupted_gene_df, known_disruptions_df
):
    filtered_df = new_serotype_no_disruptions.remove_known_disruptions(
        disrupted_gene_df, known_disruptions_df
    )
    unique_samples = new_serotype_no_disruptions.find_unique_samples(filtered_df)
    assert unique_samples == list()
