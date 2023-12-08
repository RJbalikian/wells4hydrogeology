"""
Test basic functionality before commits, deployments, etc. 

Currently a placeholder.
"""

import w4h

def test_run():
    """"Test basic functionality"""
    try:
        resource_dict = w4h.get_resources()
        print(resource_dict)
        test_passed=True
    except:
        test_passed=False
    assert test_passed
