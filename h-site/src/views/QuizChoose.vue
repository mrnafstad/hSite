<template>
	<div>
		<h1>Select a category!</h1>
		<Quiz v-if="selected != null"
			:field="fields.find(field => field.value === selected).text"
			:url="fields.find(field => field.value === selected).url" />
		<div v-else>
			<b-form-select v-model="selected" :options="fields"></b-form-select>
			<div v-if="selected != null" class="mt-3">Selected: <strong>{{ selected.text }}</strong></div>
		</div>
	</div>
</template>

<script>
import Quiz from '../components/Quiz.vue'
	export default {
		name: 'QuizChoice',
		components: {
			Quiz
		},
		data() {
			return {
				fields: [
					{
						value: null,
						text: 'Please select an option'
					},
					{
						value: 'a',
						text: 'General knowlege',
						url: 'https://opentdb.com/api.php?amount=10&type=multiple'
					},
					{
						value: 'b',
						text: 'Mythology',
						url: 'https://opentdb.com/api.php?amount=10&category=20&type=multiple'
					},
					{
						value: 'c',
						text: 'Entertainment: Books',
						url: 'https://opentdb.com/api.php?amount=10&category=10&type=multiple'
					},
					{
						value: 'd',
						text: 'Entertainment: Film',
						url: 'https://opentdb.com/api.php?amount=10&category=11&type=multiple'
					},
					{
						value: 'e',
						text: 'Entertainment: Music',
						url: 'https://opentdb.com/api.php?amount=10&category=12&type=multiple'
					},
					{
						value: 'f',
						text: 'Science and nature',
						url: 'https://opentdb.com/api.php?amount=10&category=17&type=multiple'
					},
					{
						value: 'g',
						text: 'Science: Computers',
						url: 'https://opentdb.com/api.php?amount=10&category=18&type=multiple'
					}
				],
				selected: null
			}
		}
	}
</script>