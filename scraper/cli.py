import pyfiglet
import yaml
from pathlib import Path
from config import get_supported_organisms, get_supported_features


def verify_input(user_input: str, supported_names) -> list:

    if user_input.strip() == '*':
        return supported_names

    parts = user_input.strip().split(' ')
    for part in parts:
        if part not in supported_names:
            return list()

    return parts


def print_banner() -> None:
    ascii_banner = pyfiglet.figlet_format("Ensembl scraper")
    print(ascii_banner)


def get_organisms_from_user() -> list:
    while True:

        print("Please select organisms you are interested in. Supported organisms are: "
              + str(get_supported_organisms()))
        organisms = input("Write names of organisms separated by space. Use '*' for all: ")
        organisms = verify_input(organisms, get_supported_organisms())
        if organisms:
            break
        print("Unrecognized organism. Please try again.")
        print('')

    print('')
    print("Organisms selected: " + str(organisms))

    return organisms


def get_features_from_user(organisms: list) -> dict:
    print("Now you can select feature classes for each organism.")

    user_input = dict()
    for organism in organisms:

        while True:

            print('=========================')
            print(organism)
            print('=========================')
            print("Supported feature classes are: " + str(get_supported_features(organism)))
            features = input("Write names of feature_classes separated by space. Use '*' for all: ")
            features = verify_input(features, get_supported_features(organism))
            if features:
                break
            print("Unrecognized feature class. Please try again.")
            print('')

        user_input[organism] = features

    return user_input


def get_root_dir_from_user() -> str:
    print('One more thing.')
    while True:
        root_dir = input('Please enter path to directory, where this program can create folders and store data: ')
        if Path(root_dir).exists():
            break
        print("Path doesn't exist. Please try again.")
        print('')
    return root_dir


def print_selected_settings(user_input: dict) -> None:
    print('')
    print('=========================')
    print('Selected settings')
    print('=========================')
    print(yaml.dump(user_input))
    print('')


def do_final_confirm() -> None:
    print('')
    input('All done. Confirm by any input to start.')
    print('')


def get_user_inputs():

    print_banner()
    organisms = get_organisms_from_user()
    user_input = get_features_from_user(organisms)
    print_selected_settings(user_input)
    root_dir = get_root_dir_from_user()
    do_final_confirm()

    return user_input, root_dir
