import json

def get_scan_params(sp_filename: str) -> dict:
    """Loads the scan parameters from a JSON file.

    The following is an example scan parameter file, with annotation. JSON does not support comments, so anything after '#' should be removed.

    {
        "name": "02041411.001",           # Required
        "growth": "Polished12072018",
        "sample": "D1",
        "afm_spring_constant": 39,        # Required [N/m]
        "afm_tip": "tip19",
        "thickness": 160E-9,              # Required [m]
        "ignored": false,
        "cantilevers": [                  # Required (at least one cantilever definition)
            {
                "name": "D1a2",           # Required
                "width": 2.7E-6,          # Required [m]
                "lin_ignore": 0,
                "fixed_edge": 24,         # Required (guess in pixels at fixed end location relative to left side of scan)
                "start": [26, 13],        # Required (top left corner [x, y] in pixels of cantilever)
                "end": [38, 22]           # Required (bottom right corner [x, y] in pixels of cantilever)
            },
            {
                "name": "D1a1",
                "width": 2.7E-6,
                "lin_ignore": 0,
                "fixed_edge": 26,
                "start": [28, 38],
                "end": [44, 49]
            }
        ]
    }

    """
    with open(sp_filename) as file:
        # Get the scan parameters.
        sc_params = json.load(file)

        # Some samples only have one cantilever, but we still want sc_params["cantilevers"] to be a list
        if not isinstance(sc_params["cantilevers"], list):
            sc_params["cantilevers"] = [sc_params["cantilevers"]]

    return sc_params
