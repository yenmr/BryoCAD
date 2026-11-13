"""Test keyword adaptation without installing model/GPU dependencies."""
import unittest
from peft_forward_compat import adapt_forward_for_peft

class CompatibilityTests(unittest.TestCase):
    def test_peft_optional_none_arguments(self):
        class Backbone:
            def forward(self,input_ids=None,labels=None,return_dict=None):
                return input_ids,labels,return_dict
        model=adapt_forward_for_peft(Backbone())
        self.assertEqual(model.forward(input_ids=[1,2],attention_mask=None,
            inputs_embeds=None,output_attentions=None,output_hidden_states=None,
            return_dict=True,labels=1),([1,2],1,True))

    def test_does_not_silently_discard_real_mask(self):
        class Backbone:
            def forward(self,input_ids=None):return input_ids
        model=adapt_forward_for_peft(Backbone())
        with self.assertRaisesRegex(TypeError,'attention_mask'):
            model.forward(input_ids=[1],attention_mask=[0])

    def test_preserves_supported_arguments(self):
        class Backbone:
            def forward(self,input_ids,inputs_embeds=None):return input_ids,inputs_embeds
        model=adapt_forward_for_peft(Backbone())
        self.assertEqual(model.forward([1],inputs_embeds=[2],attention_mask=None),([1],[2]))

    def test_existing_kwargs_forward_is_unchanged(self):
        class Backbone:
            def forward(self,**kwargs):return kwargs
        model=Backbone();original=model.forward
        self.assertIs(adapt_forward_for_peft(model),model)
        self.assertEqual(model.forward,original)

if __name__=='__main__':unittest.main()
