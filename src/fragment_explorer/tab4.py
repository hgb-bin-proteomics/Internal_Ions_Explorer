import streamlit as st

from .util.constants import DIV_COLOR


def main(argv=None) -> None:

    readme = ""

    with open("README.md", "r") as f:
        readme = f.read()
        f.close()

    _documentation = st.markdown(readme)

    _divider = st.subheader("", divider = DIV_COLOR)

    col_l, col_m, col_r = st.columns(3)

    with col_m:
        _ext_docs = st.link_button("Read full documentation!",
                                  url = "https://internal-ions.vercel.app/",
                                  type = "primary",
                                  width = "stretch")
