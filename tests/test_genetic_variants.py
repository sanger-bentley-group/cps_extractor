import pytest

from lib.genetic_variants import GeneticVariants


@pytest.fixture
def variant():
    genetic_variant = GeneticVariants("tests/test_data/18B.gff")
    return genetic_variant


def test_get_features(variant):
    features = variant.get_features()
    print(features)
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
        {"name": "unidentified", "product": "Glucosyl-1-phosphate transferase"},
        {"name": "unidentified", "product": "Cps23fF"},
        {"name": "unidentified", "product": "Glycosyltransferase"},
        {"name": "unidentified", "product": "Family 2 glycosyl transferase"},
        {"name": "unidentified", "product": "Glycosyl transferase"},
        {"name": "unidentified", "product": "Flippase Wzx"},
        {"name": "unidentified", "product": "O-Antigen ligase"},
        {"name": "hypothetical protein", "product": "hypothetical protein"},
        {"name": "unidentified", "product": "acyltransferase"},
        {"name": "wchX", "product": "Glycerol phosphotransferase WchX"},
        {"name": "tnp", "product": "transposase"},
        {"name": "unidentified", "product": "Heme-based aerotactic transducer"},
        {"name": "rfbA", "product": "glucose-1-phosphate thymidylyltransferase RfbA"},
        {"name": "unidentified", "product": "dTDP-4-dehydrorhamnose 3,5-epimerase"},
        {"name": "rfbB", "product": "dTDP-glucose 4,6-dehydratase"},
        {"name": "rfbD", "product": "dTDP-4-dehydrorhamnose reductase"},
        {"name": "aliA", "product": "UDP-galactopyranose mutase"},
        {"name": "glf", "product": "UDP-galactopyranose glf mutase"},
    ]
