# -*- coding: utf-8 -*-
from collections import OrderedDict


PUBLIC_ALT_SCHEME_ALIASES = OrderedDict(
    [
        ("baseline", "baseline"),
        ("fc_only_alpha2", "fc_only_alpha2"),
        ("k2_only_alpha2", "k2_only_alpha2"),
        ("fc_perc_k2_alpha6", "fc_perc_k2_alpha6"),
        ("fc_perc_k2_alpha8", "fc_perc_k2_alpha8"),
        ("fc_perc_k2_alpha10", "fc_perc_k2_alpha10"),
        ("fc_perc_k2_alpha12", "fc_perc_k2_alpha12"),
        ("fc_perc_k2_alpha14", "fc_perc_k2_alpha14"),
        ("fc_perc_k2_alpha10_clipped", "fc_perc_k2_alpha10_clipped"),
        ("physically_inverted_fc_perc_k2_alpha6", "physically_inverted_fc_perc_k2_alpha6"),
        ("physically_inverted_fc_perc_k2_alpha10", "physically_inverted_fc_perc_k2_alpha10"),
        ("forward_fc_perc_k2_alpha3", "forward_fc_perc_k2_alpha3"),
        ("forward_fc_perc_k2_exponent2", "forward_fc_perc_k2_exponent2"),
    ]
)


def public_to_legacy_scheme_names(names):
    selected = []
    for name in names:
        name = name.strip()
        if not name:
            continue
        if name in PUBLIC_ALT_SCHEME_ALIASES:
            selected.append(PUBLIC_ALT_SCHEME_ALIASES[name])
            continue
        if name in PUBLIC_ALT_SCHEME_ALIASES.values():
            selected.append(name)
            continue
        raise KeyError(
            f"Unknown ALT scheme '{name}'. Public names: {', '.join(PUBLIC_ALT_SCHEME_ALIASES.keys())}"
        )
    return selected


def legacy_to_public_scheme_name(name):
    for public_name, legacy_name in PUBLIC_ALT_SCHEME_ALIASES.items():
        if name == legacy_name:
            return public_name
    return name
