from fragment_explorer.util import converter


def test_streamlit():
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file("../streamlit_app.py", default_timeout=30.0)
    at.run(timeout=60.0)
    assert not at.exception


def test_fragment_parse():
    c = converter.JSONConverter
    assert c._parse_fragment_code('c:zdot@2:4(+1)[]') == (2, 4, 'c', 'zdot', 1, '')
    assert c._parse_fragment_code('cdot:z@2:4(+1)[H2O]') == (2, 4, 'cdot', 'z', 1, 'H2O')
    assert c._parse_fragment_code('cdot:z@2:4(+1)') == (2, 4, 'cdot', 'z', 1, '')
