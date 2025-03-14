"""
Test basic functionality before commits, deployments, etc.

Runs w4h.run() on local sample data.
"""
import w4h


def test_run():
    """"Test basic functionality"""
    try:
        resource_dict = w4h.get_resources()

        w4h.run(well_data=resource_dict['well_data'],
                model_grid=resource_dict['model_grid'],
                surf_elev_grid=resource_dict['surf_elev'],
                bedrock_elev_grid=resource_dict['bedrock_elev'],
                lith_dict=resource_dict['LithologyDict_Exact'],
                target_dict=resource_dict['LithInterps_FineCoarse'],
                study_area=resource_dict['study_area'],
                verbose=True)

        test_passed = True
    except Exception as e:
        print('ERROR:', e)
        test_passed = False

    assert test_passed


if __name__ == "__main__":
    test_run()
