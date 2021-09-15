class Item(dict):
    """
    Impletment of dict which support key adding in attribute way
    e.g. instance_dict.attr1 = 'value'
    """

    def __init__(self, *args, **kwargs):
        super(Item, self).__init__(*args, **kwargs)

    def __getattr__(self, name):
        """overriding attribute getter and setter

        Args:
            name (str): name of attribute

        Raises:
            AttributeError: Attribute not found

        Returns:
            Object: value of key
        """
        try:
            return self[name]
        except:
            raise AttributeError(
                'Item object has no attribute name {}'.format(name))
        # return super().__getattribute__(name)

    def __setattr__(self, name, value):
        """set attribute

        Args:
            name (str): name of attribute
            value (Object): value set to the attribute
        """
        self[name] = value
