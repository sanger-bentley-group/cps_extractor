import pytest
import pandas as pd

from lib.genetic_variants import GeneticVariants


@pytest.fixture
def variant():
    genetic_variant = GeneticVariants()
    return genetic_variant


def test_get_features(variant):
    features = variant.get_features("tests/test_data/18B.gff")

    assert features == [
        {
            "name": "cps2a",
            "product": "Anionic cell wall polymer biosynthesis enzyme TagV/TagU, LytR-Cps2A-Psr (LCP) family (peptidoglycan teichoic acid transferase)",
        },
        {
            "name": "cps4B",
            "product": "capsular polysaccharide biosynthesis protein Cps4B",
        },
        {
            "name": "cpsC",
            "product": "capsular polysaccharide biosynthesis protein CpsC",
        },
        {
            "name": "cps4D",
            "product": "capsular polysaccharide biosynthesis protein Cps4D",
        },
        {"name": "tnp", "product": "Glucosyl-1-phosphate transferase"},
        {"name": "tnp", "product": "Cps23fF"},
        {"name": "tnp", "product": "Glycosyltransferase"},
        {"name": "tnp", "product": "Family 2 glycosyl transferase"},
        {"name": "tnp", "product": "Glycosyl transferase"},
        {"name": "tnp", "product": "Flippase Wzx"},
        {"name": "tnp", "product": "O-Antigen ligase"},
        {"name": "tnp", "product": "hypothetical protein"},
        {"name": "tnp", "product": "acyltransferase"},
        {"name": "wchX", "product": "Glycerol phosphotransferase WchX"},
        {"name": "tnp", "product": "transposase"},
        {"name": "tnp", "product": "Heme-based aerotactic transducer"},
        {"name": "rfbA", "product": "glucose-1-phosphate thymidylyltransferase RfbA"},
        {"name": "tnp", "product": "dTDP-4-dehydrorhamnose 3,5-epimerase"},
        {"name": "rfbB", "product": "dTDP-glucose 4,6-dehydratase"},
        {"name": "rfbD", "product": "dTDP-4-dehydrorhamnose reductase"},
        {"name": "aliA", "product": "UDP-galactopyranose mutase"},
        {"name": "glf", "product": "UDP-galactopyranose glf mutase"},
    ]


def test_assign_groups(variant):
    data = [
        {"sample": "ERR311103", "genetic_group": "group1"},
        {"sample": "ERR311104", "genetic_group": "group1"},
        {"sample": "ERR311105", "genetic_group": "group1"},
        {"sample": "ERR311105", "genetic_group": "group2"},
        {"sample": "ERR311106", "genetic_group": "group2"},
        {"sample": "ERR311107", "genetic_group": "group1"},
        {"sample": "ERR311108", "genetic_group": "group3"},
    ]

    actual = pd.DataFrame(data)

    test_df = variant.assign_groups("tests/test_data/md5.csv")

    pd.testing.assert_frame_equal(actual, test_df)
