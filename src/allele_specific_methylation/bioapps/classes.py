import os
import requests


URL = os.environ.get("BIOAPPS_API_URL")
USERNAME = os.environ.get("BIOAPPS_USERNAME")
PASSWORD = os.environ.get("BIOAPPS_PASSWORD")


class RequestHandler:
    def __init__(
        self,
        url: str = URL,
        username: str = USERNAME,
        password: str = PASSWORD,
        headers: dict | None = None,
    ):
        if not url:
            msg = "BIOAPPS_API_URL environment variable is not set. Please set it to the BioApps API URL."
            raise ValueError(msg)
        if not username:
            msg = "BIOAPPS_USERNAME environment variable is not set. Please set it to your BioApps username."
            raise ValueError(msg)
        if not password:
            msg = "BIOAPPS_PASSWORD environment variable is not set. Please set it to your BioApps password."
            raise ValueError(msg)

        self.__request_handler = requests.Session()
        self.__url = url
        self.__headers = (
            {"Content-Type": "application/json", "Accept": "application/json"} if headers is None else headers
        )

        auth_request = {
            "username": username,
            "password": password,
        }
        response = self.__request_handler.post(
            os.path.join(self.__url, "session"),
            json=auth_request,
            headers=self.__headers
        )

        if response.status_code == requests.codes.ok:
            token = response.json().get("token")
            self.__headers.update({"X-Token": token})
        else:
            raise Exception(f"Authentication failed: {response.status_code} {response.text}")

    @property
    def handler(self):
        return self.__request_handler

    @property
    def headers(self):
        return self.__headers

    @property
    def url(self):
        return self.__url

    def get(self, endpoint: str, params: dict | None = None):
        response = self.__request_handler.get(
            os.path.join(self.__url, endpoint),
            params=params,
            headers=self.__headers,
        )
        if response.status_code == requests.codes.ok:
            return response
        else:
            raise Exception(f"Request failed: ({response.status_code}) {response.text}")
